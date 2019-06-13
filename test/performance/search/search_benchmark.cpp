// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <memory>
#include <random>
#include <utility>
#include <vector>


#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/algorithm/all.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <stdlib.h>

using namespace seqan3;
using namespace seqan3::test;

template <Alphabet alphabet_t>
void mutate_insertion(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha(0, alphabet_size<alphabet_t> - 1);
    std::uniform_int_distribution<size_t> random_pos(0, std::ranges::size(seq) - overlap);
    alphabet_t cbase;
    seq.insert(seq.begin() + random_pos(gen), alphabet_t{}.assign_rank(dis_alpha(gen)));
}

template <Alphabet alphabet_t>
void mutate_deletion(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> random_pos(0, std::ranges::size(seq) - overlap);
    seq.erase(seq.begin() + random_pos(gen));
}

template <Alphabet alphabet_t>
void mutate_substitution(std::vector<alphabet_t> & seq, size_t const overlap, size_t const seed)
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

template <typename alphabet_t>
auto generate_repeating_sequence(size_t const template_length = 5000,
                                 size_t const repeats = 20,
                                 double const repeats_length = 1,
                                 size_t const seed = 0)
{
    std::vector<alphabet_t> seqTemplate = generate_sequence<alphabet_t>(template_length, 0, seed);

    size_t len = template_length * repeats_length;
    int simulated_errors = 5;
    len = (len + simulated_errors  > template_length) ? template_length - simulated_errors : len;
    std::vector<std::vector<alphabet_t> > collection = generate_reads(seqTemplate, repeats, len, simulated_errors, 0.15, 0.15);

    std::vector<alphabet_t> ref;
    for(std::vector<alphabet_t> seq : collection){
        ref.insert(ref.end(), seq.begin(), seq.end());
    }

    return ref;
}

struct Options
{
    size_t const sequence_length;
    bool const repeating;
    size_t const number_of_reads;
    size_t const read_length;
    double const prob_insertion;
    double const prob_deletion;
    double const simulated_errors;
    uint8_t const searched_errors;
    uint8_t const strata;
    int repeats = 20;
};

//============================================================================
//  undirectional; trivial_search, collection, dna4, all-mapping
//============================================================================

void unidirectional_search_all_collection(benchmark::State & state, Options && o)
{
    size_t set_size = 10; 
    std::vector<std::vector<seqan3::dna4> > collection;
    std::vector<std::vector<seqan3::dna4> > reads;
    for (unsigned i = 0; i < set_size; ++i)
    {
        std::vector<seqan3::dna4> seq = generate_sequence<seqan3::dna4>(o.sequence_length, 0, 0);
        collection.push_back(seq);
        std::vector<std::vector<seqan3::dna4> > seq_reads = generate_reads(seq, o.number_of_reads, o.read_length, o.simulated_errors, o.prob_insertion, o.prob_deletion, true);
        reads.insert(reads.end(), seq_reads.begin(), seq_reads.end());
    }

    fm_index index{collection};
    configuration cfg = search_cfg::max_error{search_cfg::total{o.searched_errors}};

    for (auto _ : state)
        auto results = search(index, reads, cfg);
}


//============================================================================
//  undirectional; trivial_search, single, dna4, all-mapping
//============================================================================

void unidirectional_search_all(benchmark::State & state, Options && o)
{
    std::vector<seqan3::dna4> ref = (o.repeating) ? generate_repeating_sequence<seqan3::dna4>(o.sequence_length/o.repeats, o.repeats, 1, 0) :
                                                  generate_sequence<seqan3::dna4>(o.sequence_length, 0, 0);

    fm_index index{ref};
    std::vector<std::vector<seqan3::dna4> > reads = generate_reads(ref, o.number_of_reads, o.read_length, o.simulated_errors, o.prob_insertion, o.prob_deletion, true);
    configuration cfg = search_cfg::max_error{search_cfg::total{o.searched_errors}};

    for (auto _ : state)
        auto results = search(index, reads, cfg);
}

//============================================================================
//  bidirectional; trivial_search, single, dna4, all-mapping
//============================================================================

void bidirectional_search_all(benchmark::State & state, Options && o)
{
    std::vector<seqan3::dna4> ref = (o.repeating) ? generate_repeating_sequence<seqan3::dna4>(o.sequence_length/o.repeats, o.repeats, 1, 0) :
                                                  generate_sequence<seqan3::dna4>(o.sequence_length, 0, 0);

    bi_fm_index index{ref};
    std::vector<std::vector<seqan3::dna4> > reads = generate_reads(ref, o.number_of_reads, o.read_length, o.simulated_errors, o.prob_insertion, o.prob_deletion, true);
    configuration cfg = search_cfg::max_error{search_cfg::total{o.searched_errors}};

    for (auto _ : state)
        auto results = search(index, reads, cfg);
}

//============================================================================
//  undirectional; trivial_search, single, dna4, stratified-all-mapping
//============================================================================

void unidirectional_search_stratified(benchmark::State & state, Options && o)
{
    std::vector<seqan3::dna4> ref = (o.repeating) ? generate_repeating_sequence<seqan3::dna4>(o.sequence_length/o.repeats, o.repeats, 1, 0) :
                                                  generate_sequence<seqan3::dna4>(o.sequence_length, 0, 0);

    fm_index index{ref};
    std::vector<std::vector<seqan3::dna4> > reads = generate_reads(ref, o.number_of_reads, o.read_length, o.simulated_errors, o.prob_insertion, o.prob_deletion, true);
    configuration cfg = search_cfg::max_error{search_cfg::total{o.searched_errors}} |
                        search_cfg::mode{search_cfg::strata{o.strata}};

    for (auto _ : state)
        auto results = search(index, reads, cfg);
}

//============================================================================
//  bidirectional; trivial_search, single, dna4, stratified-all-mapping
//============================================================================

void bidirectional_search_stratified(benchmark::State & state, Options && o)
{
    std::vector<seqan3::dna4> ref = (o.repeating) ? generate_repeating_sequence<seqan3::dna4>(o.sequence_length/o.repeats, o.repeats, 1, 0) :
                                                  generate_sequence<seqan3::dna4>(o.sequence_length, 0, 0);

    bi_fm_index index{ref};
    std::vector<std::vector<seqan3::dna4> > reads = generate_reads(ref, o.number_of_reads, o.read_length, o.simulated_errors, o.prob_insertion, o.prob_deletion, true);
    configuration cfg = search_cfg::max_error{search_cfg::total{o.searched_errors}} |
                        search_cfg::mode{search_cfg::strata{o.strata}};

    for (auto _ : state)
        auto results = search(index, reads, cfg);
}


BENCHMARK_CAPTURE(unidirectional_search_all_collection, highErrorReadsSearch0, Options{10'000, false, 10, 50, 0.18, 0.18, 3, 0, 0});
BENCHMARK_CAPTURE(unidirectional_search_all_collection, highErrorReadsSearch1, Options{10'000, false, 10, 50, 0.18, 0.18, 3, 1, 0});
BENCHMARK_CAPTURE(unidirectional_search_all_collection, highErrorReadsSearch2, Options{10'000, false, 10, 50, 0.18, 0.18, 3, 2, 0});
BENCHMARK_CAPTURE(unidirectional_search_all_collection, highErrorReadsSearch3, Options{10'000, false, 10, 50, 0.18, 0.18, 3, 3, 0});

BENCHMARK_CAPTURE(unidirectional_search_all, lowErrorReadsSearch0, Options{100'000, false, 50, 50, 0.18, 0.18, 1, 0, 0});
BENCHMARK_CAPTURE(unidirectional_search_all, lowErrorReadsSearch1, Options{100'000, false, 50, 50, 0.18, 0.18, 1, 1, 0});
BENCHMARK_CAPTURE(unidirectional_search_all, lowErrorReadsSearch2, Options{100'000, false, 50, 50, 0.18, 0.18, 1, 2, 0});
BENCHMARK_CAPTURE(unidirectional_search_all, lowErrorReadsSearch3, Options{100'000, false, 50, 50, 0.18, 0.18, 1, 3, 0});
BENCHMARK_CAPTURE(unidirectional_search_all, highErrorReadsSearch0, Options{100'000, false, 50, 50, 0.18, 0.18, 3, 0, 0});
BENCHMARK_CAPTURE(unidirectional_search_all, highErrorReadsSearch1, Options{100'000, false, 50, 50, 0.18, 0.18, 3, 1, 0});
BENCHMARK_CAPTURE(unidirectional_search_all, highErrorReadsSearch2, Options{100'000, false, 50, 50, 0.18, 0.18, 3, 2, 0});
BENCHMARK_CAPTURE(unidirectional_search_all, highErrorReadsSearch3, Options{100'000, false, 50, 50, 0.18, 0.18, 3, 3, 0});
BENCHMARK_CAPTURE(unidirectional_search_all, highErrorReadsSearch0rep, Options{100'000, true, 50, 50, 0.18, 0.18, 3, 0, 0});
BENCHMARK_CAPTURE(unidirectional_search_all, highErrorReadsSearch1rep, Options{100'000, true, 50, 50, 0.18, 0.18, 3, 1, 0});
BENCHMARK_CAPTURE(unidirectional_search_all, highErrorReadsSearch2rep, Options{100'000, true, 50, 50, 0.18, 0.18, 3, 2, 0});
BENCHMARK_CAPTURE(unidirectional_search_all, highErrorReadsSearch3rep, Options{100'000, true, 50, 50, 0.18, 0.18, 3, 3, 0});

BENCHMARK_CAPTURE(bidirectional_search_all, lowErrorReadsSearch0, Options{100'000, false, 50, 50, 0.18, 0.18, 1, 0, 0});
BENCHMARK_CAPTURE(bidirectional_search_all, lowErrorReadsSearch1, Options{100'000, false, 50, 50, 0.18, 0.18, 1, 1, 0});
BENCHMARK_CAPTURE(bidirectional_search_all, lowErrorReadsSearch2, Options{100'000, false, 50, 50, 0.18, 0.18, 1, 2, 0});
BENCHMARK_CAPTURE(bidirectional_search_all, lowErrorReadsSearch3, Options{100'000, false, 50, 50, 0.18, 0.18, 1, 3, 0});
BENCHMARK_CAPTURE(bidirectional_search_all, highErrorReadsSearch0, Options{100'000, false, 50, 50, 0.18, 0.18, 3, 0, 0});
BENCHMARK_CAPTURE(bidirectional_search_all, highErrorReadsSearch1, Options{100'000, false, 50, 50, 0.18, 0.18, 3, 1, 0});
BENCHMARK_CAPTURE(bidirectional_search_all, highErrorReadsSearch2, Options{100'000, false, 50, 50, 0.18, 0.18, 3, 2, 0});
BENCHMARK_CAPTURE(bidirectional_search_all, highErrorReadsSearch3, Options{100'000, false, 50, 50, 0.18, 0.18, 3, 3, 0});
BENCHMARK_CAPTURE(bidirectional_search_all, highErrorReadsSearch0rep, Options{100'000, true, 50, 50, 0.18, 0.18, 3, 0, 0});
BENCHMARK_CAPTURE(bidirectional_search_all, highErrorReadsSearch1rep, Options{100'000, true, 50, 50, 0.18, 0.18, 3, 1, 0});
BENCHMARK_CAPTURE(bidirectional_search_all, highErrorReadsSearch2rep, Options{100'000, true, 50, 50, 0.18, 0.18, 3, 2, 0});
BENCHMARK_CAPTURE(bidirectional_search_all, highErrorReadsSearch3rep, Options{100'000, true, 50, 50, 0.18, 0.18, 3, 3, 0});

BENCHMARK_CAPTURE(unidirectional_search_stratified, lowErrorReadsSearch3Strata0rep, Options{100'000, true, 50, 50, 0.18, 0.18, 1, 3, 0});
BENCHMARK_CAPTURE(unidirectional_search_stratified, lowErrorReadsSearch3Strata1rep, Options{100'000, true, 50, 50, 0.18, 0.18, 1, 3, 1});
BENCHMARK_CAPTURE(unidirectional_search_stratified, lowErrorReadsSearch3Strata2rep, Options{100'000, true, 50, 50, 0.18, 0.18, 1, 3, 2});
BENCHMARK_CAPTURE(unidirectional_search_stratified, highErrorReadsSearch3Strata0rep, Options{100'000, true, 50, 50, 0.18, 0.18, 3, 3, 0});
BENCHMARK_CAPTURE(unidirectional_search_stratified, highErrorReadsSearch3Strata1rep, Options{100'000, true, 50, 50, 0.18, 0.18, 3, 3, 1});
BENCHMARK_CAPTURE(unidirectional_search_stratified, highErrorReadsSearch3Strata2rep, Options{100'000, true, 50, 50, 0.18, 0.18, 3, 3, 2});

BENCHMARK_CAPTURE(bidirectional_search_stratified, lowErrorReadsSearch3Strata0rep, Options{100'000, true, 50, 50, 0.18, 0.18, 1, 3, 0});
BENCHMARK_CAPTURE(bidirectional_search_stratified, lowErrorReadsSearch3Strata1rep, Options{100'000, true, 50, 50, 0.18, 0.18, 1, 3, 1});
BENCHMARK_CAPTURE(bidirectional_search_stratified, lowErrorReadsSearch3Strata2rep, Options{100'000, true, 50, 50, 0.18, 0.18, 1, 3, 2});
BENCHMARK_CAPTURE(bidirectional_search_stratified, highErrorReadsSearch3Strata0rep, Options{100'000, true, 50, 50, 0.18, 0.18, 3, 3, 0});
BENCHMARK_CAPTURE(bidirectional_search_stratified, highErrorReadsSearch3Strata1rep, Options{100'000, true, 50, 50, 0.18, 0.18, 3, 3, 1});
BENCHMARK_CAPTURE(bidirectional_search_stratified, highErrorReadsSearch3Strata2rep, Options{100'000, true, 50, 50, 0.18, 0.18, 3, 3, 2});

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
