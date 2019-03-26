// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>
#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/algorithm/all.hpp>

using namespace seqan3;

template <Alphabet alphabet_t>
auto generate_sequence_seqan3(size_t const len = 500,
                              size_t const variance = 0,
                              size_t const seed = 0)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha(0, alphabet_size_v<alphabet_t> - 1);
    std::uniform_int_distribution<size_t> dis_length(len - variance, len + variance);

    std::vector<alphabet_t> sequence;

    size_t length = dis_length(gen);
    for (size_t l = 0; l < length; ++l)
        sequence.push_back(alphabet_t{}.assign_rank(dis_alpha(gen)));

    return sequence;
}

template <Alphabet alphabet_t>
void mutate_insertion(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed = 0)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha(0, alphabet_size_v<alphabet_t> - 1);
    std::uniform_int_distribution<size_t> random_pos(0, std::ranges::size(seq) - overlap);
    alphabet_t cbase;
    seq.insert(seq.begin() + random_pos(gen), alphabet_t{}.assign_rank(dis_alpha(gen)));
}

template <Alphabet alphabet_t>
void mutate_deletion(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed = 0)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> random_pos(0, std::ranges::size(seq) - overlap);
    seq.erase(seq.begin() + random_pos(gen));
}

template <Alphabet alphabet_t>
void mutate_substitution(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed = 0)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha_short(0, alphabet_size_v<alphabet_t> - 2);
    std::uniform_int_distribution<size_t> random_pos(0, std::ranges::size(seq) - overlap);
    alphabet_t & cbase = seq[random_pos(gen)];
    uint8_t crank = to_rank(cbase);
    uint8_t rrank = dis_alpha_short(gen);
    if (rrank >=  crank)
        ++rrank;
    cbase.assign_rank(rrank);
}

template <Alphabet alphabet_t>
auto generate_reads(std::vector<alphabet_t> & ref,
                    size_t const number_of_reads,
                    size_t const read_length,
                    size_t const simulated_errors,
                    float const prob_insertion,
                    float const prob_deletion,
                    size_t const seed = 0)
{
    std::vector<std::vector<alphabet_t> > reads;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> seeds (0, SIZE_MAX);
    std::uniform_int_distribution<size_t> random_pos(0, std::ranges::size(ref) - read_length - simulated_errors);
    std::uniform_real_distribution<double> probability(0.0, 1.0);
    for (size_t i = 0; i < number_of_reads; ++i)
    {
        size_t rpos = random_pos(gen);
        std::vector<alphabet_t> read_tmp{ref.begin() + rpos,
            ref.begin() + rpos + read_length + simulated_errors};
        for (size_t j = 0; j < simulated_errors; ++j)
        {
            double prob = probability(gen);
            //Substitution
            if (prob_insertion + prob_deletion < prob)
            {
                mutate_substitution(read_tmp, simulated_errors, seeds(gen));
            }
            //Insertion
            else if (prob_insertion < prob)
            {
                mutate_insertion(read_tmp, simulated_errors, seeds(gen));
            }
            //Deletion
            else
            {
                mutate_deletion(read_tmp, simulated_errors, seeds(gen));
            }
        }
        read_tmp.erase(read_tmp.begin() + read_length, read_tmp.end());
        reads.push_back(read_tmp);
    }
    return reads;
}

//============================================================================
//  undirectional; trivial_search, single, dna4, searched_error = 3
//============================================================================

void unidirectional_search(benchmark::State & state)
{
    uint8_t const simulated_errors = state.range(0);
    uint8_t const searched_errors = state.range(1);
    size_t reference_length{100'000};
    int number_of_reads{100};
    size_t read_length{100};
    float prob_insertion{0.18};
    float prob_deletion{0.18};

    std::vector<seqan3::dna4> ref = generate_sequence_seqan3<seqan3::dna4>(reference_length);
    fm_index<std::vector<seqan3::dna4> > index{ref};
    std::vector<std::vector<seqan3::dna4> > reads = generate_reads(ref, number_of_reads, read_length, simulated_errors, prob_insertion, prob_deletion);
    configuration cfg = search_cfg::max_error{search_cfg::total{searched_errors}};

    for (auto _ : state)
    {
        auto results = search(index, reads, cfg);
    }
}

BENCHMARK(unidirectional_search)
    ->Args({0, 1})
    ->Args({1, 1})
    ->Args({0, 3})
    ->Args({1, 3})
    ->Args({2, 3})
    ->Args({3, 3});

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
