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
//  insert left to right
// ============================================================================
template <typename gap_decorator_t, bool gapped_flag>
void insert_left2right(benchmark::State & state)
{
    unsigned int seq_len = state.range(0);
    using size_type = typename gap_decorator_t::size_type;
    using sequence_type = remove_cvref_t<detail::unaligned_seq_t<gap_decorator_t>>;
    sequence_type seq(seq_len, 'A'_dna4);

    // vector of sampled gap lengths for each position
    std::vector<size_type> gaps(seq_len, 0);

    // determine sum of gaps and non-gap symbols for not exceeding targeted sequence length
    if (gapped_flag)
    {
        sample<size_type>(gaps, seq_len, state.range(1)/100.0);
        resize<size_type, sequence_type>(gaps, seq, seq_len);
    }
    // initialize with (truncated) sequence and insert gaps from left to right
    gap_decorator_t gap_decorator;
    assign_unaligned(gap_decorator, seq);

    // insert gaps before starting benchmark
    if (gapped_flag)
        insert_gaps<gap_decorator_t>(gaps, gap_decorator, seq_len);

    auto it = std::ranges::begin(gap_decorator);
    for (auto _ : state)
    {
        if (it == std::ranges::end(gap_decorator) || it == (--std::ranges::end(gap_decorator)))
            it = std::ranges::begin(gap_decorator);
        it = insert_gap(gap_decorator, it, 1);
        ++++it;
    }
}

// Insert gaps of length 1 from left to right into UNGAPPED sequence
BENCHMARK_TEMPLATE(insert_left2right, gap_decorator<const std::vector<dna4> &>, false)->Apply(custom_arguments);
BENCHMARK_TEMPLATE(insert_left2right, std::vector<gapped<dna4>>, false)->Apply(custom_arguments);

// Insert gaps of length 1 from left to right into GAPPED sequence
BENCHMARK_TEMPLATE(insert_left2right, gap_decorator<const std::vector<dna4> &>, true)->Apply(custom_arguments);
BENCHMARK_TEMPLATE(insert_left2right, std::vector<gapped<dna4>>, true)->Apply(custom_arguments);

// ============================================================================
//  insert right to left
// ============================================================================
template <typename gap_decorator_t, bool gapped_flag>
void insert_right2left(benchmark::State & state)
{
    unsigned int seq_len = state.range(0);
    assert(seq_len > 0);
    using size_type = typename gap_decorator_t::size_type;
    using sequence_type = remove_cvref_t<detail::unaligned_seq_t<gap_decorator_t>>;
    sequence_type seq(seq_len, 'A'_dna4);
    // vector of sample<size_type>d gap lengths for each position
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
    // end() != begin() due to assert above on sequence length
    auto it = std::ranges::end(gap_decorator);
    for (auto _ : state)
    {
        it = insert_gap(gap_decorator, it, 1);
        state.PauseTiming();
        if (it == std::ranges::begin(gap_decorator))
            it = std::ranges::end(gap_decorator);
        --it;
        state.ResumeTiming();
    }
}

// Insert gaps of length 1 from left to right into UNGAPPED sequence
BENCHMARK_TEMPLATE(insert_right2left, gap_decorator<const std::vector<dna4> &>, false)->Apply(custom_arguments);
BENCHMARK_TEMPLATE(insert_right2left, std::vector<gapped<dna4>>, false)->Apply(custom_arguments);
// Insert gaps of length 1 from left to right into GAPPED sequence
BENCHMARK_TEMPLATE(insert_right2left, gap_decorator<const std::vector<dna4> &>, true)->Apply(custom_arguments);
BENCHMARK_TEMPLATE(insert_right2left, std::vector<gapped<dna4>>, true)->Apply(custom_arguments);

BENCHMARK_MAIN();
