// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
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

// ============================================================================
//  insert at random position
// ============================================================================
template <typename gap_decorator_t, bool gapped_flag>
    requires aligned_sequence<gap_decorator_t>
void insert_random(benchmark::State & state)
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
    std::vector<size_type> access_positions(1 << 10);
    std::generate(access_positions.begin(), access_positions.end(),
        [&](){return uni_dis(generator);});
    size_t j = 0;
    for (auto _ : state)
    {
        state.PauseTiming();
        auto it = std::ranges::begin(gap_decorator);
        std::ranges::advance(it, access_positions[j++]);
        state.ResumeTiming();

        insert_gap(gap_decorator, it, 1);
        j %= 1 << 10;
    }
}

// Insert gaps of length 1 at random position into UNGAPPED sequence
BENCHMARK_TEMPLATE(insert_random, gap_decorator<const std::vector<dna4> &>, false)->Apply(custom_arguments);
BENCHMARK_TEMPLATE(insert_random, std::vector<gapped<dna4>>, false)->Apply(custom_arguments);
// Insert gaps of length 1 at random position into GAPPED sequence
BENCHMARK_TEMPLATE(insert_random, gap_decorator<const std::vector<dna4> &>, true)->Apply(custom_arguments);
BENCHMARK_TEMPLATE(insert_random, std::vector<gapped<dna4>>, true)->Apply(custom_arguments);

// ============================================================================
//  delete at random position
// ============================================================================
template <typename gap_decorator_t, bool gapped_flag>
    requires aligned_sequence<gap_decorator_t>
void delete_random(benchmark::State & state)
{
    unsigned int seq_len = state.range(0);
    using iterator_type = std::ranges::iterator_t<gap_decorator_t>;
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

    // presample delete positions
    std::vector<size_t> access_positions(1 << 10);
    std::generate(access_positions.begin(), access_positions.end(),
        [&](){return uni_dis(generator);});

    iterator_type first, last;
    size_t j = 0;
    for (auto _ : state)
    {
        state.PauseTiming();
        first = std::ranges::begin(gap_decorator);
        std::ranges::advance(first, access_positions[j++]);
        first = insert_gap(gap_decorator, first, 2);
        last = first;
        std::ranges::advance(last, 2);
        j %= 1 << 10;
        state.ResumeTiming();
        erase_gap(gap_decorator, first, last);
    }
}

// Erase gaps at random position from initially GAPPED sequence
BENCHMARK_TEMPLATE(delete_random, gap_decorator<const std::vector<dna4> &>, true)->Apply(custom_arguments);
BENCHMARK_TEMPLATE(delete_random, std::vector<gapped<dna4>>, true)->Apply(custom_arguments);

BENCHMARK_MAIN();
