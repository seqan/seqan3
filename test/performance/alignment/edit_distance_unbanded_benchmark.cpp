// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <iostream>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/alignment/align_pairwise_edit_distance.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/utility/views/zip.hpp>

#ifdef SEQAN3_HAS_SEQAN2
#    include <seqan/align.h>
#    include <seqan/basic.h>
#    include <seqan/find.h>
#    include <seqan/sequence.h>
#endif

constexpr auto edit_distance_cfg =
    seqan3::align_cfg::method_global{} | seqan3::align_cfg::edit_scheme | seqan3::align_cfg::output_score{};

// Shortcut to determine the alignment result type.
template <typename seq1_t, typename seq2_t, typename config_t>
using alignment_result_type_t =
    seqan3::alignment_result<typename seqan3::detail::align_result_selector<seq1_t, seq2_t, config_t>::type>;

// ============================================================================
//  edit_distance; score; dna4; single
// ============================================================================

void seqan3_edit_distance_dna4(benchmark::State & state)
{
    using seqan3::test::edit_distance_algorithm;

    size_t sequence_length = 500;
    auto seq1 = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 0);
    auto seq2 = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 1);
    int score = 0;

    auto algorithm = edit_distance_algorithm::select<decltype(seq1), decltype(seq2)>(edit_distance_cfg);

    for (auto _ : state)
        score += algorithm(seq1, seq2, edit_distance_cfg).score();

    state.counters["score"] = score;
    state.counters["cells"] =
        seqan3::test::pairwise_cell_updates(std::views::single(std::tie(seq1, seq2)), edit_distance_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

void seqan3_edit_distance_dna4_selector(benchmark::State & state)
{
    size_t sequence_length = 500;
    auto seq1 = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 0);
    auto seq2 = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 1);
    int score = 0;

    for (auto _ : state)
    {
        for (auto && rng : align_pairwise(std::tie(seq1, seq2), edit_distance_cfg))
            score += rng.score();
    }

    state.counters["score"] = score;
    state.counters["cells"] =
        seqan3::test::pairwise_cell_updates(std::views::single(std::tie(seq1, seq2)), edit_distance_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

#ifdef SEQAN3_HAS_SEQAN2
void seqan2_edit_distance_dna4(benchmark::State & state)
{
    using seqan3::test::edit_distance_algorithm_seqan2;

    size_t sequence_length = 500;
    auto seq1 = seqan3::test::generate_sequence_seqan2<seqan2::Dna>(sequence_length, 0, 0);
    auto seq2 = seqan3::test::generate_sequence_seqan2<seqan2::Dna>(sequence_length, 0, 1);
    int score = 0;

    seqan3::configuration align_cfg = seqan3::align_cfg::method_global{};
    auto algorithm_seqan2 = edit_distance_algorithm_seqan2::select<decltype(seq1), decltype(seq2)>(align_cfg);

    for (auto _ : state)
        score += algorithm_seqan2(seq1, seq2);

    state.counters["score"] = score;
    state.counters["cells"] =
        seqan3::test::pairwise_cell_updates(std::views::single(std::tie(seq1, seq2)), edit_distance_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

void seqan2_edit_distance_generic_dna4(benchmark::State & state)
{
    size_t sequence_length = 500;
    auto seq1 = seqan3::test::generate_sequence_seqan2<seqan2::Dna>(sequence_length, 0, 0);
    auto seq2 = seqan3::test::generate_sequence_seqan2<seqan2::Dna>(sequence_length, 0, 1);
    int score = 0;

    for (auto _ : state)
        score += seqan2::globalAlignmentScore(seq1, seq2, seqan2::Score<int>{0, -1, -1});

    state.counters["score"] = score;
    state.counters["cells"] =
        seqan3::test::pairwise_cell_updates(std::views::single(std::tie(seq1, seq2)), edit_distance_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  edit_distance; score; dna4; set
// ============================================================================

void seqan3_edit_distance_dna4_collection(benchmark::State & state)
{
    using seqan3::test::edit_distance_algorithm;

    size_t sequence_length = 500;
    size_t set_size = 100;

    auto vec = seqan3::test::generate_sequence_pairs<seqan3::dna4>(sequence_length, set_size);
    int score = 0;

    auto algorithm = edit_distance_algorithm::select<decltype(std::get<0>(vec[0])), decltype(std::get<1>(vec[0]))>(
        edit_distance_cfg);

    for (auto _ : state)
    {
        for (auto && [seq1, seq2] : vec)
            score += algorithm(seq1, seq2, edit_distance_cfg).score();
    }

    state.counters["score"] = score;
    state.counters["cells"] = seqan3::test::pairwise_cell_updates(vec, edit_distance_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

void seqan3_edit_distance_dna4_collection_selector(benchmark::State & state)
{
    size_t sequence_length = 500;
    size_t set_size = 100;

    auto vec = seqan3::test::generate_sequence_pairs<seqan3::dna4>(sequence_length, set_size);
    int score = 0;

    for (auto _ : state)
    {
        for (auto && rng : align_pairwise(vec, edit_distance_cfg))
            score += rng.score();
    }

    state.counters["score"] = score;
    state.counters["cells"] = seqan3::test::pairwise_cell_updates(vec, edit_distance_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

#ifdef SEQAN3_HAS_SEQAN2
void seqan2_edit_distance_dna4_collection(benchmark::State & state)
{
    using seqan3::test::edit_distance_algorithm_seqan2;

    size_t sequence_length = 500;
    size_t set_size = 100;

    auto [vec1, vec2] = seqan3::test::generate_sequence_pairs_seqan2<seqan2::Dna>(sequence_length, set_size);
    int score = 0;

    seqan3::configuration align_cfg = seqan3::align_cfg::method_global{};
    auto algorithm_seqan2 = edit_distance_algorithm_seqan2::select<decltype(vec1[0]), decltype(vec2[0])>(align_cfg);

    for (auto _ : state)
    {
        for (unsigned i = 0; i < set_size; ++i)
            score += algorithm_seqan2(vec1[i], vec2[i]);
    }

    state.counters["score"] = score;
    state.counters["cells"] = seqan3::test::pairwise_cell_updates(seqan3::views::zip(vec1, vec2), edit_distance_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

void seqan2_edit_distance_dna4_generic_collection(benchmark::State & state)
{
    size_t sequence_length = 500;
    size_t set_size = 100;

    auto [vec1, vec2] = seqan3::test::generate_sequence_pairs_seqan2<seqan2::Dna>(sequence_length, set_size);
    int score = 0;

    for (auto _ : state)
    {
        for (int score_ : seqan2::globalAlignmentScore(vec1, vec2, seqan2::Score<int>{0, -1, -1}))
            score += score_;
    }

    state.counters["score"] = score;
    state.counters["cells"] = seqan3::test::pairwise_cell_updates(seqan3::views::zip(vec1, vec2), edit_distance_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK(seqan3_edit_distance_dna4);
BENCHMARK(seqan3_edit_distance_dna4_selector);
#ifdef SEQAN3_HAS_SEQAN2
BENCHMARK(seqan2_edit_distance_dna4);
BENCHMARK(seqan2_edit_distance_generic_dna4);
#endif
BENCHMARK(seqan3_edit_distance_dna4_collection);
BENCHMARK(seqan3_edit_distance_dna4_collection_selector);
#ifdef SEQAN3_HAS_SEQAN2
BENCHMARK(seqan2_edit_distance_dna4_collection);
BENCHMARK(seqan2_edit_distance_dna4_generic_collection);
#endif

BENCHMARK_MAIN();
