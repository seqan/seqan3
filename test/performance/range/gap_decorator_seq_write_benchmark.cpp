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

#include "gap_decorator_helper.hpp"

using seqan3::operator""_dna4;

using gap_sequence_gap_decorator = seqan3::gap_decorator<std::vector<seqan3::dna4> const &>;
using gap_sequence_vector = std::vector<seqan3::gapped<seqan3::dna4>>;

// ============================================================================
//  insert left to right
// ============================================================================
template <typename gap_decorator_t, bool gapped_flag>
void insert_left2right(benchmark::State & state)
{
    unsigned int seq_len = state.range(0);
    using size_type = typename gap_decorator_t::size_type;
    using sequence_type = std::remove_cvref_t<seqan3::detail::unaligned_seq_t<gap_decorator_t>>;
    sequence_type seq(seq_len, 'A'_dna4);

    // vector of sampled gap lengths for each position
    std::vector<size_type> gaps(seq_len, 0);

    // determine sum of gaps and non-gap symbols for not exceeding targeted sequence length
    if (gapped_flag)
    {
        sample<size_type>(gaps, seq_len, state.range(1) / 100.0);
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
        ++ ++it;
    }
}

// Insert gaps of length 1 from left to right into UNGAPPED sequence
BENCHMARK_TEMPLATE(insert_left2right, gap_sequence_gap_decorator, false)->Apply(custom_arguments);
BENCHMARK_TEMPLATE(insert_left2right, gap_sequence_vector, false)->Apply(custom_arguments);

// Insert gaps of length 1 from left to right into GAPPED sequence
BENCHMARK_TEMPLATE(insert_left2right, gap_sequence_gap_decorator, true)->Apply(custom_arguments);
BENCHMARK_TEMPLATE(insert_left2right, gap_sequence_vector, true)->Apply(custom_arguments);

// ============================================================================
//  insert right to left
// ============================================================================
template <typename gap_decorator_t, bool gapped_flag>
void insert_right2left(benchmark::State & state)
{
    unsigned int seq_len = state.range(0);
    assert(seq_len > 0);
    using size_type = typename gap_decorator_t::size_type;
    using sequence_type = std::remove_cvref_t<seqan3::detail::unaligned_seq_t<gap_decorator_t>>;
    sequence_type seq(seq_len, 'A'_dna4);
    // vector of sample<size_type>d gap lengths for each position
    std::vector<size_type> gaps(seq_len, 0);

    // determine sum of gaps and non-gap symbols for not exceeding targeted sequence length
    if constexpr (gapped_flag)
    {
        sample<size_type>(gaps, seq_len, state.range(1) / 100.0);
        resize<size_type, sequence_type>(gaps, seq, seq_len);
    }
    // initialize with (truncated) sequence and insert gaps from left to right
    gap_decorator_t gap_decorator;
    assign_unaligned(gap_decorator, seq);

    // insert gaps before starting benchmark
    if constexpr (gapped_flag)
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
BENCHMARK_TEMPLATE(insert_right2left, gap_sequence_gap_decorator, false)->Apply(custom_arguments);
BENCHMARK_TEMPLATE(insert_right2left, gap_sequence_vector, false)->Apply(custom_arguments);
// Insert gaps of length 1 from left to right into GAPPED sequence
BENCHMARK_TEMPLATE(insert_right2left, gap_sequence_gap_decorator, true)->Apply(custom_arguments);
BENCHMARK_TEMPLATE(insert_right2left, gap_sequence_vector, true)->Apply(custom_arguments);

BENCHMARK_MAIN();
