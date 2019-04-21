// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

// written for the purpose of the Bch Softwarepraktikum, early 2019 FU Berlin
// @author: Clemens Cords <clemenscords@fu-berlin.de>

#include <benchmark/benchmark.h>
#include <cstring>
#include <random>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

using namespace seqan3;

std::size_t const LENGTH = 50;
unsigned int const SEED = 1234;

// ============================================================================
// construct as many indices as possible
// ============================================================================
template<typename index_t>
void construct_index(benchmark::State & state)
{
    using alphabet_t = typename index_t::char_type;

    // generate sequences in advance to not skew benchmarking results
    std::vector<std::vector<alphabet_t>> sequences;
    std::size_t n_sequences = 10000; // arbitrary

    for (std::size_t i = 0; i < n_sequences; ++i)
        sequences.push_back(test::generate_sequence<alphabet_t>(LENGTH, 0, SEED));

    std::size_t i = 0;

    for (auto _ : state)
    {
        i = i < n_sequences + 1 ? i + 1 : 0;
        index_t index{sequences[i]};
    }
}

// TODO <Clemens C.>: seqan2 comparison

BENCHMARK_TEMPLATE(construct_index, fm_index<std::vector<dna4>>);
BENCHMARK_TEMPLATE(construct_index, fm_index<std::vector<aa27>>);
BENCHMARK_TEMPLATE(construct_index, fm_index<std::vector<phred63>>);
BENCHMARK_TEMPLATE(construct_index, fm_index<std::vector<char>>);

BENCHMARK_TEMPLATE(construct_index, bi_fm_index<std::vector<dna4>>);
BENCHMARK_TEMPLATE(construct_index, bi_fm_index<std::vector<aa27>>);
BENCHMARK_TEMPLATE(construct_index, bi_fm_index<std::vector<phred63>>);
BENCHMARK_TEMPLATE(construct_index, bi_fm_index<std::vector<char>>);

BENCHMARK_MAIN();
