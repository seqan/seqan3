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

#include "gap_decorator_helper.hpp"

using namespace seqan3;

// Apply benchmarks with custom ranges for grid parameters sequence length and gap proportions
static void CustomArguments(benchmark::internal::Benchmark* b) {
    std::array<long long int, 5> gap_percentages = {1, 5, 25, 50, 75};
    for (long long int seq_len = 4; seq_len <= (1 << 14); seq_len <<= 2)
    {
        for (auto gap_percentage : gap_percentages)
            b->Args({seq_len, gap_percentage});
    }
}

// ============================================================================
//  insert at random position
// ============================================================================
template <typename gap_decorator_t, bool gapped_flag>
static void insert_random(benchmark::State& state)
{
    unsigned int seq_len = state.range(0);
    using size_type = typename gap_decorator_t::size_type;
    using sequence_type = remove_cvref_t<detail::unaligned_seq_t<gap_decorator_t>>;
    sequence_type seq(seq_len, 'A'_dna4);

    // vector of sampled gap lengths for each position
    std::vector<size_type> gaps(seq_len, 0);

    // determine sum of gaps and non-gap symbols for not exceeding targeted sequence length
    if constexpr(gapped_flag)
    {
        sample<size_type>(gaps, seq_len, state.range(1)/100.0);
        resize<size_type, sequence_type>(gaps, seq, seq_len);
    }

    // initialize with (truncated) sequence and insert gaps from left to right
    gap_decorator_t gap_decorator;
    assign_unaligned(gap_decorator, seq);

    // insert gaps before starting benchmark
    if constexpr(gapped_flag)
        insert_gaps<gap_decorator_t>(gaps, gap_decorator, seq_len);

    std::mt19937 generator(time(0)); //Standard mersenne_twister_engine seeded with current time
    std::uniform_real_distribution<> uni_dis{0.0, static_cast<double>(seq_len)};

    // presample insert positions
    std::vector<size_t> access_positions;
    access_positions.resize(1 << 18);
    for (size_t i = 0; i < access_positions.size(); ++i)
        access_positions[i] = uni_dis(generator);
    size_t j = 0;
    for (auto _ : state)
    {
        state.PauseTiming();
        auto it = std::ranges::begin(gap_decorator);
        std::ranges::advance(it, access_positions[j++]);
        if (j == (1 << 18))
            j = 0;
        state.ResumeTiming();
        insert_gap(gap_decorator, it, 1);
    }
}

// Insert gaps of length 1 at random position into UNGAPPED sequence
BENCHMARK_TEMPLATE(insert_random, gap_decorator_anchor_set<const std::vector<dna4> &>, false)->Apply(CustomArguments);
BENCHMARK_TEMPLATE(insert_random, std::vector<gapped<dna4>>, false)->Apply(CustomArguments);
// Insert gaps of length 1 at random position into GAPPED sequence
BENCHMARK_TEMPLATE(insert_random, gap_decorator_anchor_set<const std::vector<dna4> &>, true)->Apply(CustomArguments);
BENCHMARK_TEMPLATE(insert_random, std::vector<gapped<dna4>>, true)->Apply(CustomArguments);

// ============================================================================
//  delete at random position
// ============================================================================
template <typename gap_decorator_t, bool gapped_flag>
static void delete_random(benchmark::State& state)
{
    unsigned int seq_len = state.range(0);
    using size_type = typename gap_decorator_t::size_type;
    using sequence_type = remove_cvref_t<detail::unaligned_seq_t<gap_decorator_t>>;
    sequence_type seq(seq_len, 'A'_dna4);

    // vector of sampled gap lengths for each position
    std::vector<size_type> gaps(seq_len, 0);

    // determine sum of gaps and non-gap symbols for not exceeding targeted sequence length
    if constexpr(gapped_flag)
    {
        sample<size_type>(gaps, seq_len, state.range(1)/100.0);
        resize<size_type, sequence_type>(gaps, seq, seq_len);
    }

    // initialize with (truncated) sequence and insert gaps from left to right
    gap_decorator_t gap_decorator;
    assign_unaligned(gap_decorator, seq);

    // insert gaps before starting benchmark
    if constexpr(gapped_flag)
        insert_gaps<gap_decorator_t>(gaps, gap_decorator, seq_len);

    std::mt19937 generator(time(0)); //Standard mersenne_twister_engine seeded with current time
    std::uniform_real_distribution<> uni_dis{0.0, static_cast<double>(seq_len)};

    std::vector<size_t> access_positions;
    access_positions.resize(1 << 18);
    for (size_t i = 0; i < access_positions.size(); ++i)
        access_positions[i] = uni_dis(generator);
    size_t j = 0;
    for (auto _ : state)
    {
        state.PauseTiming();
        auto first = std::ranges::begin(gap_decorator);
        std::ranges::advance(first, access_positions[j++]);
        first = insert_gap(gap_decorator, first, 2);
        auto last = std::ranges::begin(gap_decorator);
        last = first;
        std::ranges::advance(last, 2);
        if (j == (1 << 18))
            j = 0;
        state.ResumeTiming();
        erase_gap(gap_decorator, first, last);
    }
}

// Erase gaps at random position from initially GAPPED sequence
BENCHMARK_TEMPLATE(delete_random, gap_decorator_anchor_set<const std::vector<dna4> &>, true)->Apply(CustomArguments);
BENCHMARK_TEMPLATE(delete_random, std::vector<gapped<dna4>>, true)->Apply(CustomArguments);

BENCHMARK_MAIN();
