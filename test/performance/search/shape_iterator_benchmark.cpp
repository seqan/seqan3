// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------


#include <random>
#include <vector>
#include <algorithm>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/search/kmer_index/shape_iterator.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

using namespace seqan3;

// performance test on dna4 type alphabet only
static void shape_iterator_hashing(benchmark::State& state)
{
    auto seq = seqan3::test::generate_sequence<seqan3::dna4>(state.range(0), 0, 0);
    shape s(8);
    auto it_end(seq.end());
    for (auto _ : state)
        for (shape_iterator it(seq.begin(), s); it != it_end; ++it)
        { }
}

static void shape_iterator_hashing_gapped(benchmark::State& state)
{
    auto seq = seqan3::test::generate_sequence<seqan3::dna4>(state.range(0), 0, 0);
    shape s(1,1,1,1,0,1,1,1);
    auto it_end(seq.end());
    for (auto _ : state)
        for (shape_iterator it(seq.begin(), s); it != it_end; ++it)
        { }
}

static void shape_iterator_hashing_random_access(benchmark::State& state)
{
    auto seq = seqan3::test::generate_sequence<seqan3::dna4>(state.range(0), 0, 0);
    shape s(8);
    auto it_end(seq.end());
    shape_iterator it(seq.begin(), s);
    for (auto _ : state)
    {
        for (int i{0}; i < state.range(0); ++i)
        {
            it[i];
        }
    }
}

static void shape_iterator_hashing_random_access_gapped(benchmark::State& state)
{
    auto seq = seqan3::test::generate_sequence<seqan3::dna4>(state.range(0), 0, 0);
    shape s(1,1,1,1,0,1,1,1);
    auto it_end(seq.end());
    shape_iterator it(seq.begin(), s);
    for (auto _ : state)
    {
        for (int i{0}; i < state.range(0); ++i)
        {
            it[i];
        }
    }
}

BENCHMARK(shape_iterator_hashing)->Range(8, 8<<12);
BENCHMARK(shape_iterator_hashing_gapped)->Range(8, 8<<12);
BENCHMARK(shape_iterator_hashing_random_access)->Range(8, 8<<12);
BENCHMARK(shape_iterator_hashing_random_access_gapped)->Range(8, 8<<12);


BENCHMARK_MAIN();
