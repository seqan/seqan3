// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <chrono>
#include <cmath>
#include <cstring>
#include <random>
#include <utility>

#include <benchmark/benchmark.h>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alignment/exception.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/range/decorator/all.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

typedef long double time_type;

/* Helper function to sample gap length for each ungapped sequence position.
 *
 * Parameter:
 * size_type    size type of underlying sequence range
 *
 * Arguments:
 * gap_vector   reference to empty vector storing the sampled gap lengths
 * size         final (aligned) sequence size
 */
template<typename size_type>
void sample(std::vector<size_type> & gap_vector, size_type size)
{
    std::default_random_engine generator;
    generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_real_distribution<double> uni(0.0, 1.0);
    std::array<double,10> cumsum = {0.6395, 0.8263, 0.8871, 0.9257, 0.9544, 0.9709, 0.9813, 0.9890, 0.9955, 1.0000};
    for (size_type i = 0; i < size; ++i){
        double y = uni(generator);
        gap_vector[i] = y;
        auto it = std::find_if(cumsum.begin(), cumsum.end(), [y](double i){return y <= i;});
        gap_vector[i] = it - cumsum.begin();
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
template<typename size_type, typename sequence_type>
void resize(std::vector<size_type> & gaps, sequence_type & seq, unsigned int seq_len)
{
    size_type letter_acc = 0;
    size_type gap_pos = 0;
    size_type gap_acc = 0;

    while (gap_pos < gaps.size() && gap_acc + letter_acc < seq_len)
    {
        if (!gaps[gap_pos])
            ++letter_acc;
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
                gap_acc += gaps[gap_pos];
        }
        ++gap_pos;
    }
    seq.resize(std::max<size_type>(1, letter_acc));  // resize ungapped sequence
    gaps.resize(gap_pos);      // trim sampled gap vector
}

/* Helper function to prepare a gapped sequence for the benchmark (case gap_flag=true)
 *
 * Parameters:
 * size_type        size type of the sequence
 * gap_decorator_t  gap decorator type, e.g. gap_decorator_anchor_set
 *
 * Arguments:
 * gaps             reference to gap vector
 * gap_decorator    reference to gap decorator
 */
template<typename size_type, typename gap_decorator_t>
void insert_gaps(std::vector<size_type> & gaps, gap_decorator_t & gap_decorator)
{
    [[maybe_unused]] size_type gap_acc = 0;
    for (size_type i = 0; i < gaps.size(); ++i)
    {
        if (gaps[i])
        {
            auto it = std::ranges::begin(gap_decorator);
            std::ranges::advance(it, std::min(i + gap_acc, gap_decorator.size()));
            insert_gap(gap_decorator, it, gaps[i]);
        }
        gap_acc += gaps[i];
    }
}

// ============================================================================
//  read left to right (looped in case #ops exceeds sequence length)
// ============================================================================
/* Parameters:
 * gap_decorator_t      gap decorator class, e.g. gap_decorator_anchor_set
 * gapped_flag          operate on already gapped (true) or ungapped sequence (false)
 */
template <typename gap_decorator_t, bool gapped_flag>
static void read_left2right(benchmark::State& state)
{
    // get target sequence length from current range state
    unsigned int seq_len = state.range(0);
    using size_type = typename gap_decorator_t::size_type;
    using sequence_type = remove_cvref_t<detail::unaligned_seq_t<gap_decorator_t>>;
    sequence_type seq(seq_len, 'A'_dna4);

    // vector of sampled gap lengths for each position
    std::vector<size_type> gaps(seq_len, 0);
    sample<size_type>(gaps, seq_len);

    // determine sum of gaps and non-gap symbols for not exceeding targeted sequence length
    if constexpr(gapped_flag)
        resize<size_type, sequence_type>(gaps, seq, seq_len);

    // initialize with (truncated) sequence and insert gaps from left to right
    gap_decorator_t gap_decorator;
    assign_unaligned(gap_decorator, seq);

    // insert gaps before starting benchmark
    if constexpr(gapped_flag)
        insert_gaps<size_type, gap_decorator_t>(gaps, gap_decorator);

    size_t op_ctr = 0;
    for (auto _ : state)
    {
        state.PauseTiming();
        size_t pos = op_ctr % seq_len;
        state.ResumeTiming();
        benchmark::DoNotOptimize(gap_decorator[pos]);
        state.PauseTiming();
        ++op_ctr;
        state.ResumeTiming();
    }
    state.counters["read_op"] = op_ctr;
}

// 1 a) Read from left to right in ungapped sequence
BENCHMARK_TEMPLATE(read_left2right, gap_decorator_anchor_set<const std::vector<dna4> &>, false)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(read_left2right, std::vector<gapped<dna4>>, false)->Range(1<<2, 1<<15);
// 1 b) Read from left to right in gapped sequence
BENCHMARK_TEMPLATE(read_left2right, gap_decorator_anchor_set<const std::vector<dna4> &>, true)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(read_left2right, std::vector<gapped<dna4>>, true)->Range(1<<2, 1<<15);

// ============================================================================
//  insert left to right
// ============================================================================
template <typename gap_decorator_t, bool gapped_flag>
static void insert_left2right(benchmark::State& state)
{
    unsigned int seq_len = state.range(0);
    using size_type = typename gap_decorator_t::size_type;
    using sequence_type = remove_cvref_t<detail::unaligned_seq_t<gap_decorator_t>>;
    sequence_type seq(seq_len, 'A'_dna4);

    // vector of sampled gap lengths for each position
    std::vector<size_type> gaps(seq_len, 0);
    sample<size_type>(gaps, seq_len);

    // determine sum of gaps and non-gap symbols for not exceeding targeted sequence length
    if constexpr(gapped_flag)
        resize<size_type, sequence_type>(gaps, seq, seq_len);

    // initialize with (truncated) sequence and insert gaps from left to right
    gap_decorator_t gap_decorator;
    assign_unaligned(gap_decorator, seq);

    // insert gaps before starting benchmark
    if constexpr(gapped_flag)
        insert_gaps<size_type, gap_decorator_t>(gaps, gap_decorator);

    size_t op_ctr = 0;
    for (auto _ : state)
    {
        state.PauseTiming();
        size_t pos = op_ctr % seq_len;
        auto it = std::ranges::begin(gap_decorator);
        std::ranges::advance(it, pos);
        state.ResumeTiming();
        insert_gap(gap_decorator, it, 1);
        state.PauseTiming();
        ++op_ctr;
        state.ResumeTiming();
    }
    state.counters["insert_op"] = op_ctr;
}

// 2 a) Insert gaps of length 1 from left to right into ungapped sequence
BENCHMARK_TEMPLATE(insert_left2right, gap_decorator_anchor_set<const std::vector<dna4> &>, false)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(insert_left2right, std::vector<gapped<dna4>>, false)->Range(1<<2, 1<<15);
// 2 b) Insert gaps of length 1 from left to right into gapped sequence
BENCHMARK_TEMPLATE(insert_left2right, gap_decorator_anchor_set<const std::vector<dna4> &>, true)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(insert_left2right, std::vector<gapped<dna4>>, true)->Range(1<<2, 1<<15);

// ============================================================================
//  insert right to left
// ============================================================================
template <typename gap_decorator_t, bool gapped_flag>
static void insert_right2left(benchmark::State& state)
{
    unsigned int seq_len = state.range(0);
    using size_type = typename gap_decorator_t::size_type;
    using sequence_type = remove_cvref_t<detail::unaligned_seq_t<gap_decorator_t>>;
    sequence_type seq(seq_len, 'A'_dna4);

    // vector of sampled gap lengths for each position
    std::vector<size_type> gaps(seq_len, 0);
    sample<size_type>(gaps, seq_len);

    // determine sum of gaps and non-gap symbols for not exceeding targeted sequence length
    if constexpr(gapped_flag)
        resize<size_type, sequence_type>(gaps, seq, seq_len);

    // initialize with (truncated) sequence and insert gaps from left to right
    gap_decorator_t gap_decorator;
    assign_unaligned(gap_decorator, seq);

    // insert gaps before starting benchmark
    if constexpr(gapped_flag)
        insert_gaps<size_type, gap_decorator_t>(gaps, gap_decorator);

    size_t op_ctr = 0;
    for (auto _ : state)
    {
        state.PauseTiming();
        size_t pos = seq_len - (op_ctr % seq_len) - 1;
        auto it = std::ranges::begin(gap_decorator);
        std::ranges::advance(it, pos);
        state.ResumeTiming();
        insert_gap(gap_decorator, it, 1);
        state.PauseTiming();
        ++op_ctr;
        state.ResumeTiming();
    }
    state.counters["insert_op"] = op_ctr;
}

// 3 a) Insert gaps of length 1 from left to right into ungapped sequence
BENCHMARK_TEMPLATE(insert_right2left, gap_decorator_anchor_set<const std::vector<dna4> &>, false)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(insert_right2left, std::vector<gapped<dna4>>, false)->Range(1<<2, 1<<15);
// 3 b) Insert gaps of length 1 from left to right into gapped sequence
BENCHMARK_TEMPLATE(insert_right2left, gap_decorator_anchor_set<const std::vector<dna4> &>, true)->Range(1<<2, 1<<15);
BENCHMARK_TEMPLATE(insert_right2left, std::vector<gapped<dna4>>, true)->Range(1<<2, 1<<15);

BENCHMARK_MAIN();
