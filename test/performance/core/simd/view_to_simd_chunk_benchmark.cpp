// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <deque>
#include <iterator>
#include <list>
#include <vector>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/view_to_simd.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>
#include <seqan3/test/performance/sequence_generator.hpp>

// ============================================================================
//  naive implementation without condition inside of hot loop
// ============================================================================

template <typename container_t, typename simd_t>
void to_simd_naive_wo_condition(benchmark::State& state)
{
    constexpr size_t simd_length = seqan3::simd::simd_traits<simd_t>::length;
    // Preparing the sequences
    std::vector<container_t> sequences;
    sequences.resize(simd_length);

    for (size_t i = 0; i < simd_length; ++i)
        std::ranges::copy(seqan3::test::generate_sequence<seqan3::dna4>(500, 10),
                          std::ranges::back_inserter(sequences[i]));

    size_t value = 0;
    for (auto _ : state)
    {
        // First sort the sequences by their lengths, but only use a proxy.
        auto sorted_sequences =
            std::views::transform(seqan3::views::zip(sequences, std::views::iota(0u, simd_length)), [] (auto && tpl)
            {
                return std::pair{std::ranges::size(std::get<0>(tpl)), std::get<1>(tpl)};
            })
            | seqan3::views::to<std::vector<std::pair<size_t, size_t>>>;

        std::ranges::sort(sorted_sequences);

        // Prepare the simd representation and transform the set.
        std::vector<simd_t, seqan3::aligned_allocator<simd_t, sizeof(simd_t)>> v;
        v.resize(sorted_sequences.back().first, seqan3::fill<simd_t>(seqan3::alphabet_size<seqan3::dna4>));

        size_t start_offset = 0;
        for (size_t k = 0; k < sorted_sequences.size(); ++k)
        {
            for (size_t i = start_offset; i < sorted_sequences[k].first; ++i)
                for (size_t j = k; j < simd_length; ++j)
                    v[i][sorted_sequences[j].second] = seqan3::to_rank(sequences[sorted_sequences[j].second][i]);

            start_offset = sorted_sequences[k].first;
        }

        for (simd_t & vec : v)
            value += vec[0];
    }

    state.counters["value"] = value;
}

// runs with contiguous_range
BENCHMARK_TEMPLATE(to_simd_naive_wo_condition, std::vector<seqan3::dna4>, seqan3::simd::simd_type_t<int8_t>);
BENCHMARK_TEMPLATE(to_simd_naive_wo_condition, std::vector<seqan3::dna4>, seqan3::simd::simd_type_t<int16_t>);
BENCHMARK_TEMPLATE(to_simd_naive_wo_condition, std::vector<seqan3::dna4>, seqan3::simd::simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(to_simd_naive_wo_condition, std::vector<seqan3::dna4>, seqan3::simd::simd_type_t<int64_t>);

BENCHMARK_TEMPLATE(to_simd_naive_wo_condition, std::deque<seqan3::dna4>, seqan3::simd::simd_type_t<int8_t>);
BENCHMARK_TEMPLATE(to_simd_naive_wo_condition, std::deque<seqan3::dna4>, seqan3::simd::simd_type_t<int16_t>);
BENCHMARK_TEMPLATE(to_simd_naive_wo_condition, std::deque<seqan3::dna4>, seqan3::simd::simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(to_simd_naive_wo_condition, std::deque<seqan3::dna4>, seqan3::simd::simd_type_t<int64_t>);

// ============================================================================
//  naive implementation with condition inside of hot loop
// ============================================================================

template <typename container_t, typename simd_t>
void to_simd_naive_w_condition(benchmark::State& state)
{
    constexpr size_t simd_length = seqan3::simd::simd_traits<simd_t>::length;
    // Preparing the sequences
    std::vector<container_t> sequences;
    sequences.resize(simd_length);

    for (size_t i = 0; i < simd_length; ++i)
        std::ranges::copy(seqan3::test::generate_sequence<seqan3::dna4>(500, 10),
                          std::ranges::back_inserter(sequences[i]));

    size_t value = 0;
    for (auto _ : state)
    {
        size_t max_size = std::ranges::size(*(std::ranges::max_element(sequences, [] (auto const & lhs,
                                                                                      auto const & rhs)
        {
            return std::ranges::size(lhs) < std::ranges::size(rhs);
        })));

        std::vector<simd_t, seqan3::aligned_allocator<simd_t, sizeof(simd_t)>> v;
        v.resize(max_size, seqan3::fill<simd_t>(seqan3::alphabet_size<seqan3::dna4>));

        for (size_t i = 0; i < max_size; ++i)
            for (size_t j = 0; j < simd_length; ++j)
                v[i][j] = (i < sequences[j].size()) ? seqan3::to_rank(sequences[j][i]) : 0;

        for (simd_t & vec : v)
            value += vec[0];
    }

    state.counters["value"] = value;
}

// runs with contiguous_range
BENCHMARK_TEMPLATE(to_simd_naive_w_condition, std::vector<seqan3::dna4>, seqan3::simd::simd_type_t<int8_t>);
BENCHMARK_TEMPLATE(to_simd_naive_w_condition, std::vector<seqan3::dna4>, seqan3::simd::simd_type_t<int16_t>);
BENCHMARK_TEMPLATE(to_simd_naive_w_condition, std::vector<seqan3::dna4>, seqan3::simd::simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(to_simd_naive_w_condition, std::vector<seqan3::dna4>, seqan3::simd::simd_type_t<int64_t>);

BENCHMARK_TEMPLATE(to_simd_naive_w_condition, std::deque<seqan3::dna4>, seqan3::simd::simd_type_t<int8_t>);
BENCHMARK_TEMPLATE(to_simd_naive_w_condition, std::deque<seqan3::dna4>, seqan3::simd::simd_type_t<int16_t>);
BENCHMARK_TEMPLATE(to_simd_naive_w_condition, std::deque<seqan3::dna4>, seqan3::simd::simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(to_simd_naive_w_condition, std::deque<seqan3::dna4>, seqan3::simd::simd_type_t<int64_t>);

// ============================================================================
//  view implementation
// ============================================================================

template <typename container_t, typename simd_t>
void to_simd(benchmark::State& state)
{
    // Preparing the sequences
    std::vector<container_t> sequences;
    sequences.resize(seqan3::simd::simd_traits<simd_t>::length);

    for (size_t i = 0; i < seqan3::simd::simd_traits<simd_t>::length; ++i)
        std::ranges::copy(seqan3::test::generate_sequence<seqan3::dna4>(500, 10),
                          std::ranges::back_inserter(sequences[i]));

    size_t value = 0;
    for (auto _ : state)
    {
        for (auto && chunk : sequences | seqan3::views::to_simd<simd_t>)
            for (simd_t const & vec : chunk)
                value += vec[0];
    }

    state.counters["value"] = value;
}

// runs with contiguous_range
BENCHMARK_TEMPLATE(to_simd, std::vector<seqan3::dna4>, seqan3::simd::simd_type_t<int8_t>);
BENCHMARK_TEMPLATE(to_simd, std::vector<seqan3::dna4>, seqan3::simd::simd_type_t<int16_t>);
BENCHMARK_TEMPLATE(to_simd, std::vector<seqan3::dna4>, seqan3::simd::simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(to_simd, std::vector<seqan3::dna4>, seqan3::simd::simd_type_t<int64_t>);

BENCHMARK_TEMPLATE(to_simd, std::deque<seqan3::dna4>, seqan3::simd::simd_type_t<int8_t>);
BENCHMARK_TEMPLATE(to_simd, std::deque<seqan3::dna4>, seqan3::simd::simd_type_t<int16_t>);
BENCHMARK_TEMPLATE(to_simd, std::deque<seqan3::dna4>, seqan3::simd::simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(to_simd, std::deque<seqan3::dna4>, seqan3::simd::simd_type_t<int64_t>);

// runs without contiguous_range
BENCHMARK_TEMPLATE(to_simd, std::list<seqan3::dna4>, seqan3::simd::simd_type_t<int8_t>);
BENCHMARK_TEMPLATE(to_simd, std::list<seqan3::dna4>, seqan3::simd::simd_type_t<int16_t>);
BENCHMARK_TEMPLATE(to_simd, std::list<seqan3::dna4>, seqan3::simd::simd_type_t<int32_t>);
BENCHMARK_TEMPLATE(to_simd, std::list<seqan3::dna4>, seqan3::simd::simd_type_t<int64_t>);

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
