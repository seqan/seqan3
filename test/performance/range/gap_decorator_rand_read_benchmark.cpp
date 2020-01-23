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
//  read at random position
// ============================================================================
template <typename gap_decorator_t, bool gapped_flag>
void read_random(benchmark::State & state)
{
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

    std::mt19937 generator(time(0)); //Standard mersenne_twister_engine seeded with current time
    std::uniform_real_distribution<> uni_dis{0.0, static_cast<double>(seq_len)};

    // sample read positions in advance
    std::vector<size_t> access_positions(1 << 10);
    std::generate(access_positions.begin(), access_positions.end(),
        [&](){return uni_dis(generator);});

    size_t j = 0, k;
    for (auto _ : state)
    {
        for (k = 0; k < 10; ++k)
            benchmark::DoNotOptimize(gap_decorator[access_positions[j + k]]);
        ++j;
        j %= (1 << 10) - 10;
    }
}

using gap_sequence_gap_decorator = seqan3::gap_decorator<const std::vector<seqan3::dna4> &>;
using gap_sequence_vector = std::vector<seqan3::gapped<seqan3::dna4>>;

// Read at random position in UNGAPPED sequence
BENCHMARK_TEMPLATE(read_random, gap_sequence_gap_decorator, false)->Apply(custom_arguments);
BENCHMARK_TEMPLATE(read_random, gap_sequence_vector, false)->Apply(custom_arguments);
// Read at random position in GAPPED sequence
BENCHMARK_TEMPLATE(read_random, gap_sequence_gap_decorator, true)->Apply(custom_arguments);
BENCHMARK_TEMPLATE(read_random, gap_sequence_vector, true)->Apply(custom_arguments);

BENCHMARK_MAIN();
