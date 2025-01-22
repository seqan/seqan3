// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <chrono>
#include <cmath>
#include <cstring>
#include <random>
#include <ranges>
#include <utility>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alignment/decorator/gap_decorator.hpp>
#include <seqan3/alignment/exception.hpp>
#include <seqan3/alphabet/all.hpp>

#define SEQAN3_LEN_LONG 1 << 18
#define SEQAN3_LEN_SHORT 1 << 12

/*
 * Apply benchmarks with custom ranges for grid parameters sequence length and gap
 * proportions. When setting the SEQAN3_LONG_TESTS macro, the longer benchmark
 * version will be executed.
 *
 *  | flag                | sequence lengths | gap proportions
 *  | --------------------|------------------|---------------------
 *  | !SEQAN3_LONG_TESTS  |  1 << [4:2:12]   | [1, 5, 50] %
 *  | SEQAN3_LONG_TESTS   |  1 << [4:2:18]   | [1, 5, 25, 50, 75] %
 */
void custom_arguments(benchmark::internal::Benchmark * b)
{
    std::vector<long long int> gap_percentages{1, 5, 50};
    long long int seq_len_max = SEQAN3_LEN_SHORT;

#ifdef SEQAN3_LONG_TESTS
    seq_len_max = SEQAN3_LEN_LONG;
    gap_percentages = {1, 5, 25, 50, 75};
#endif

    for (long long int seq_len = 16; seq_len <= seq_len_max; seq_len <<= 2)
    {
        for (auto gap_percentage : gap_percentages)
            b->Args({seq_len, gap_percentage});
    }
}

/* Helper function to sample gap length for each ungapped sequence position.
 * Distribution is derived from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC419611.
 * Parameter:
 * size_type    size type of underlying sequence range
 *
 * Arguments:
 * gap_vector   reference to empty vector storing the sampled gap lengths
 * size         final (aligned) sequence size
 * gap_density  targeted relative amount of gap symbols relative to the final size
 */
template <typename size_type>
void sample(std::vector<size_type> & gap_vector, size_type size, double gap_density)
{
    std::default_random_engine generator;
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> uni{0.0, 1.0};
    // cdf for gap length distribution, probability for length 0 is 64%, for length 1 83%-64% and so forth
    // the bin index i corresponds to the gap length and height(i)-height(i-1) to its probability
    std::array<double, 10> cumsum = {0.6395, 0.8263, 0.8871, 0.9257, 0.9544, 0.9709, 0.9813, 0.9890, 0.9955, 1.0000};
    size_type gap_acc = 0;
    for (size_type i = 0; i < size; ++i)
    {
        double y = uni(generator);
        auto it = std::find_if(cumsum.begin(),
                               cumsum.end(),
                               [y](double bar)
                               {
                                   return y <= bar;
                               });
        gap_vector[i] = it - cumsum.begin();
        gap_acc += gap_vector[i];
        if (gap_acc >= gap_density * size)
            break;
    }
}

/* Helper function to adjust the ungapped sequence length w.r.t. sampled gaps s.t.
 * the gapped sequence length does not exceed the targeted length.
 *
 * Parameters:
 * size_type        size type of the sequence
 * sequence_type    ungapped sequence type
 *
 * Arguments:
 * gaps         reference to the gap vector
 * seq          reference to ungapped sequence
 * seq_len      final sequence length
 */
template <typename size_type, typename sequence_type>
void resize(std::vector<size_type> & gaps, sequence_type & seq, unsigned int seq_len)
{
    size_type letter_acc = 0;
    size_type gap_pos = 0;
    size_type gap_acc = 0;

    while (gap_pos < gaps.size() && gap_acc + letter_acc < seq_len)
    {
        if (!gaps[gap_pos])
        {
            ++letter_acc;
        }
        else
        {
            if (letter_acc + gap_acc + gaps[gap_pos] > seq_len)
            {
                gaps[gap_pos] = seq_len - gap_acc - letter_acc;
                gap_acc += gaps[gap_pos];
                ++gap_pos;
                break;
            }
            else
            {
                gap_acc += gaps[gap_pos];
            }
        }
        ++gap_pos;
    }
    seq.resize(std::max<size_type>(1, letter_acc)); // resize ungapped sequence
    gaps.resize(gap_pos);                           // trim sampled gap vector
}

/* Helper function to prepare a gapped sequence for the benchmark (case gap_flag=true)
 *
 * Parameters:
 * gap_decorator_t  gap decorator type, e.g. gap_decorator
 *
 * Arguments:
 * gaps             reference to gap vector
 * gap_decorator    reference to gap decorator
 */
template <typename gap_decorator_t>
void insert_gaps(std::vector<typename gap_decorator_t::size_type> & gaps,
                 gap_decorator_t & gap_decorator,
                 typename gap_decorator_t::size_type target_len)
{
    typename gap_decorator_t::size_type gap_acc = 0;
    typename gap_decorator_t::size_type insert_pos = 0;
    std::ranges::iterator_t<gap_decorator_t> it;
    for (typename gap_decorator_t::size_type i = 0; i < gaps.size(); ++i)
    {
        if (gaps[i])
        {
            it = std::ranges::begin(gap_decorator);
            insert_pos = std::min(i + gap_acc, gap_decorator.size());
            std::ranges::advance(it, insert_pos);
            insert_gap(gap_decorator, it, gaps[i]);
            if (insert_pos + gaps[i] >= target_len)
                return;
        }
        gap_acc += gaps[i];
    }
}
