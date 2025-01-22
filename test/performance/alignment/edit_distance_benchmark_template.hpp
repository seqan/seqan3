// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <benchmark/benchmark.h>

#include <seqan3/test/alignment/align_pairwise_edit_distance.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/test/seqan2.hpp>

// ----------------------------------------------------------------------------
//  seqan3 edit distance pairwise alignment benchmarks
// ----------------------------------------------------------------------------

// This directly benchmarks the alignment algorithm without going through the whole align_pairwise function
template <typename sequence_pair_generator_t, typename align_configs_t>
void seqan3_align_pairwise_edit_distance_benchmark(benchmark::State & state,
                                                   sequence_pair_generator_t sequence_pair_generator,
                                                   align_configs_t && edit_distance_cfg)
{
    using seqan3::test::edit_distance_algorithm;

    auto sequence_pair_or_pairs = sequence_pair_generator(state);
    constexpr bool collection_benchmark = sequence_pair_generator.is_collection;
    using sequence_t = typename sequence_pair_generator_t::sequence_t;

    int64_t total = 0;

    auto algorithm = edit_distance_algorithm::select<sequence_t, sequence_t>(edit_distance_cfg);

    if constexpr (collection_benchmark)
    {
        for (auto _ : state)
            for (auto && [sequence1, sequence2] : sequence_pair_or_pairs)
                total += algorithm(sequence1, sequence2, edit_distance_cfg).score();
    }
    else
    {
        for (auto _ : state)
        {
            auto & [sequence1, sequence2] = sequence_pair_or_pairs;
            total += algorithm(sequence1, sequence2, edit_distance_cfg).score();
        }
    }

    std::conditional_t<collection_benchmark, decltype(std::views::all), decltype(std::views::single)> view_adaptor{};
    state.counters["cells"] =
        seqan3::test::pairwise_cell_updates(view_adaptor(sequence_pair_or_pairs), edit_distance_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

#ifdef SEQAN3_HAS_SEQAN2

// ----------------------------------------------------------------------------
//  seqan2 edit distance pairwise alignment benchmarks
// ----------------------------------------------------------------------------

template <typename sequence_pair_generator_t, typename align_cfg_t>
void seqan2_align_pairwise_edit_distance_benchmark(benchmark::State & state,
                                                   sequence_pair_generator_t sequence_pair_generator,
                                                   align_cfg_t seqan3_align_cfg)
{
    using seqan3::test::edit_distance_algorithm_seqan2;

    auto [sequences1, sequences2] = sequence_pair_generator(state);
    constexpr bool collection_benchmark = sequence_pair_generator.is_collection;
    using sequence_t = typename sequence_pair_generator_t::sequence_t;

    auto algorithm_seqan2 = edit_distance_algorithm_seqan2::select<sequence_t, sequence_t>(seqan3_align_cfg);

    int64_t total = 0;

    if constexpr (collection_benchmark)
    {
        for (auto _ : state)
            for (size_t i = 0u; i < seqan2::length(sequences1); ++i)
                total += algorithm_seqan2(sequences1[i], sequences2[i]);
    }
    else
    {
        for (auto _ : state)
            total += algorithm_seqan2(sequences1, sequences2);
    }

    std::conditional_t<collection_benchmark, decltype(std::views::all), decltype(std::views::single)> view_adaptor{};
    auto sequence_pairs_view = seqan3::views::zip(view_adaptor(sequences1), view_adaptor(sequences2));
    state.counters["cells"] = seqan3::test::pairwise_cell_updates(sequence_pairs_view, seqan3_align_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

#endif // SEQAN3_HAS_SEQAN2
