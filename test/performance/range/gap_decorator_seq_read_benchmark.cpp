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
    for (long long int seq_len = 4; seq_len <= (1 << 18); seq_len <<= 2)
    {
        for (auto gap_percentage : gap_percentages)
            b->Args({seq_len, gap_percentage});
    }
}

// ============================================================================
//  read left to right (looped in case #ops exceeds sequence length)
// ============================================================================
/* Parameters:
 * gap_decorator_t      gap decorator class, e.g. gap_decorator_anchor_set
 * gapped_flag          operate on already gapped (true) or ungapped sequence (false)
 */
template <typename gap_decorator_t, bool gapped_flag>
static void read_left2right(benchmark::State& state)
{
    // get target sequence length from current range state
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

    size_t op_ctr = 0;

    for (auto _ : state)
    {
        size_t pos = op_ctr % seq_len;
        benchmark::DoNotOptimize(gap_decorator[pos]);
        ++op_ctr;
    }
}

// Read from left to right in UNGAPPED sequence
BENCHMARK_TEMPLATE(read_left2right, gap_decorator_anchor_set<const std::vector<dna4> &>, false)->Apply(CustomArguments);
BENCHMARK_TEMPLATE(read_left2right, std::vector<gapped<dna4>>, false)->Apply(CustomArguments);
// Read from left to right in GAPPED sequence
BENCHMARK_TEMPLATE(read_left2right, gap_decorator_anchor_set<const std::vector<dna4> &>, true)->Apply(CustomArguments);
BENCHMARK_TEMPLATE(read_left2right, std::vector<gapped<dna4>>, true)->Apply(CustomArguments);

BENCHMARK_MAIN();
