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

using seqan3::operator""_dna4;

// ============================================================================
//  read left to right (looped in case #ops exceeds sequence length)
// ============================================================================
template <typename gap_decorator_t, bool gapped_flag>
void read_left2right(benchmark::State & state)
{
    // get target sequence length from current range state
    unsigned int seq_len = state.range(0);
    using size_type = typename gap_decorator_t::size_type;
    using sequence_type = seqan3::remove_cvref_t<seqan3::detail::unaligned_seq_t<gap_decorator_t>>;
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

    auto it = std::ranges::begin(gap_decorator);
    for (auto _ : state)
    {
        benchmark::DoNotOptimize(*it++);
        if (it == std::ranges::end(gap_decorator))
            it = std::ranges::begin(gap_decorator);
    }
}

using gap_sequence_gap_decorator = seqan3::gap_decorator<const std::vector<seqan3::dna4> &>;
using gap_sequence_vector = std::vector<seqan3::gapped<seqan3::dna4>>;

// Read from left to right in UNGAPPED sequence
BENCHMARK_TEMPLATE(read_left2right, gap_sequence_gap_decorator, false)->Apply(custom_arguments);
BENCHMARK_TEMPLATE(read_left2right, gap_sequence_vector, false)->Apply(custom_arguments);
// Read from left to right in GAPPED sequence
BENCHMARK_TEMPLATE(read_left2right, gap_sequence_gap_decorator, true)->Apply(custom_arguments);
BENCHMARK_TEMPLATE(read_left2right, gap_sequence_vector, true)->Apply(custom_arguments);

BENCHMARK_MAIN();
