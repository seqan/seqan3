// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <iostream>
#include <memory>
#include <random>
#include <utility>
#include <vector>


#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/algorithm/all.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/io/all.hpp>
#include <stdlib.h>




using namespace seqan3;
using namespace seqan3::test;

template <Alphabet alphabet_t>
void mutate_insertion(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed = 0)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha(0, alphabet_size<alphabet_t> - 1);
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
    std::uniform_int_distribution<uint8_t> dis_alpha_short(0, alphabet_size<alphabet_t> - 2);
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
                    float const simulated_errors_v,
                    float const prob_insertion,
                    float const prob_deletion,
                    bool const error_dis = false,
                    size_t const seed = 0)
{
    std::vector<std::vector<alphabet_t> > reads;
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> seeds (0, SIZE_MAX);
    std::normal_distribution<> dis_error_count{0, simulated_errors_v};
    std::uniform_real_distribution<double> probability(0.0, 1.0);
    for (size_t i = 0; i < number_of_reads; ++i)
    {
        //TODO avoid mutating the same position multiple times

        int simulated_errors;
        if(!error_dis)
            simulated_errors = std::round(simulated_errors_v);
        else
            simulated_errors = abs(std::round(dis_error_count(gen)));
//         debug_stream << i << ": " << simulated_errors << "\n";

        std::uniform_int_distribution<size_t> random_pos(0, std::ranges::size(ref) - read_length - simulated_errors);
        size_t rpos = random_pos(gen);
        std::vector<alphabet_t> read_tmp{ref.begin() + rpos,
            ref.begin() + rpos + read_length + simulated_errors};
        for (int j = 0; j < simulated_errors; ++j)
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
//  undirectional; trivial_search, single, dna4, all-mapping
//============================================================================

void unidirectional_search_all(benchmark::State & state)
{

    size_t const sequence_length = state.range(0);
    size_t const number_of_reads = state.range(1);
    size_t const read_length = state.range(2);
    double const prob_insertion = static_cast<double>(state.range(3))/100;
    double const prob_deletion = static_cast<double>(state.range(4))/100;
    double const simulated_errors = state.range(5);
    uint8_t const searched_errors = state.range(6);


    std::vector<seqan3::dna4> ref = generate_sequence<seqan3::dna4>(sequence_length, 0, 0);
    fm_index<std::vector<seqan3::dna4> > index{ref};
    std::vector<std::vector<seqan3::dna4> > reads = generate_reads(ref, number_of_reads, read_length, simulated_errors, prob_insertion, prob_deletion);
    configuration cfg = search_cfg::max_error{search_cfg::total{searched_errors}};
    for (auto _ : state)
    {
        auto results = search(index, reads, cfg);
    }
}

//============================================================================
//  bidirectional; trivial_search, single, dna4, all-mapping
//============================================================================

void bidirectional_search_all(benchmark::State & state)
{
    double const simulated_errors = state.range(0);
    uint8_t const searched_errors = state.range(1);
    size_t sequence_length{100'000};
    int number_of_reads{100};
    size_t read_length{100};
    float prob_insertion{0.18};
    float prob_deletion{0.18};

    std::vector<seqan3::dna4> ref = generate_sequence<seqan3::dna4>(sequence_length, 0, 0);
    bi_fm_index<std::vector<dna4> > index{ref};
    std::vector<std::vector<seqan3::dna4> > reads = generate_reads(ref, number_of_reads, read_length, simulated_errors, prob_insertion, prob_deletion);
    configuration cfg = search_cfg::max_error{search_cfg::total{searched_errors}};

    for (auto _ : state)
    {
        auto results = search(index, reads, cfg);
    }
}

//============================================================================
//  undirectional; trivial_search, single, dna4, stratified-all-mapping
//============================================================================

void unidirectional_search_stratified(benchmark::State & state)
{
    size_t const sequence_length = state.range(0);
    size_t const number_of_reads = state.range(1);
    size_t const read_length = state.range(2);
    double const prob_insertion = static_cast<double>(state.range(3))/100;
    double const prob_deletion = static_cast<double>(state.range(4))/100;
    uint8_t const simulated_errors = static_cast<double>(state.range(5))/100;
    uint8_t const searched_errors = state.range(6);
    uint8_t const strata = state.range(7);

    std::vector<seqan3::dna4> ref = generate_sequence<seqan3::dna4>(sequence_length, 0, 0);
    fm_index<std::vector<seqan3::dna4> > index{ref};
    std::vector<std::vector<seqan3::dna4> > reads = generate_reads(ref, number_of_reads, read_length, simulated_errors, prob_insertion, prob_deletion, true);
    configuration cfg = search_cfg::max_error{search_cfg::total{searched_errors}} |
                        search_cfg::mode{search_cfg::strata{strata}};

    for (auto _ : state)
    {
        auto results = search(index, reads, cfg);
    }
}

void bidirectional_search_stratified(benchmark::State & state)
{
    size_t const sequence_length = state.range(0);
    size_t const number_of_reads = state.range(1);
    size_t const read_length = state.range(2);
    double const prob_insertion = static_cast<double>(state.range(3))/100;
    double const prob_deletion = static_cast<double>(state.range(4))/100;
    uint8_t const simulated_errors = static_cast<double>(state.range(5))/100;
    uint8_t const searched_errors = state.range(6);
    uint8_t const strata = state.range(7);

    std::vector<seqan3::dna4> ref = generate_sequence<seqan3::dna4>(sequence_length, 0, 0);
    bi_fm_index<std::vector<seqan3::dna4> > index{ref};
    std::vector<std::vector<seqan3::dna4> > reads = generate_reads(ref, number_of_reads, read_length, simulated_errors, prob_insertion, prob_deletion, true);
    configuration cfg = search_cfg::max_error{search_cfg::total{searched_errors}} |
                        search_cfg::mode{search_cfg::strata{strata}};

    for (auto _ : state)
    {
        auto results = search(index, reads, cfg);
    }
}


BENCHMARK(unidirectional_search_all)
    ->Args({100'000, 100, 100, 18, 18, 0, 1})
    ->Args({100'000, 100, 100, 18, 18, 1, 1})
    ->Args({100'000, 100, 100, 18, 18, 0, 3})
    ->Args({100'000, 100, 100, 18, 18, 1, 3})
    ->Args({100'000, 100, 100, 18, 18, 2, 3})
    ->Args({100'000, 100, 100, 18, 18, 3, 3})
    ->Args({100'000, 100, 100, 18, 18, 3, 3});

BENCHMARK(bidirectional_search_all)
    ->Args({0, 1})
    ->Args({1, 1})
    ->Args({0, 3})
    ->Args({1, 3})
    ->Args({2, 3})
    ->Args({3, 3});


BENCHMARK(unidirectional_search_stratified)
    ->Args({100'000, 100, 100, 18, 18, 100, 3, 0})
    ->Args({100'000, 100, 100, 18, 18, 100, 3, 1})
    ->Args({100'000, 100, 100, 18, 18, 100, 3, 2})
    ->Args({100'000, 100, 100, 18, 18, 200, 3, 0})
    ->Args({100'000, 100, 100, 18, 18, 200, 3, 1})
    ->Args({100'000, 100, 100, 18, 18, 200, 3, 2});

BENCHMARK(bidirectional_search_stratified)
    ->Args({100'000, 100, 100, 18, 18, 100, 3, 0})
    ->Args({100'000, 100, 100, 18, 18, 100, 3, 1})
    ->Args({100'000, 100, 100, 18, 18, 100, 3, 2})
    ->Args({100'000, 100, 100, 18, 18, 200, 3, 0})
    ->Args({100'000, 100, 100, 18, 18, 200, 3, 1})
    ->Args({100'000, 100, 100, 18, 18, 200, 3, 2});

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
