// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <iostream>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include <range/v3/view/zip.hpp>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/seqan2.hpp>

#ifdef SEQAN3_HAS_SEQAN2
#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/find.h>
#include <seqan/sequence.h>
#endif

using namespace seqan3;
using namespace seqan3::test;

constexpr auto edit_distance_cfg = align_cfg::edit |
                                   align_cfg::result{with_score};

#ifdef SEQAN3_HAS_SEQAN2
template <typename alphabet_t>
int global_edit_distance_seqan2(
    [[maybe_unused]] seqan::String<alphabet_t> & text,
    [[maybe_unused]] seqan::String<alphabet_t> & needle)
{
    using text_t = seqan::String<alphabet_t>;
    using needle_t = text_t;
    using pattern_t = seqan::Pattern<needle_t, seqan::MyersUkkonenGlobal>;

    pattern_t pattern(needle, std::numeric_limits<int>::max());
    seqan::Finder<text_t> finder(text);
    while (seqan::find(finder, pattern));

    return -static_cast<int>(pattern.errors);
}
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  edit_distance; score; dna4; single
// ============================================================================

void seqan3_edit_distance_dna4(benchmark::State & state)
{
    size_t sequence_length = 500;
    auto seq1 = generate_sequence<seqan3::dna4>(sequence_length, 0, 0);
    auto seq2 = generate_sequence<seqan3::dna4>(sequence_length, 0, 1);
    int score = 0;

    for (auto _ : state)
    {
        detail::edit_distance_unbanded edit_distance{seq1, seq2, edit_distance_cfg};
        edit_distance(0u);
        score += edit_distance.score().value_or(0u);
    }

    state.counters["score"] = score;
    state.counters["cells"] = pairwise_cell_updates(ranges::view::single(std::tie(seq1, seq2)), edit_distance_cfg);
    state.counters["CUPS"] = cell_updates_per_second(state.counters["cells"]);
}

void seqan3_edit_distance_dna4_selector(benchmark::State & state)
{
    size_t sequence_length = 500;
    auto seq1 = generate_sequence<seqan3::dna4>(sequence_length, 0, 0);
    auto seq2 = generate_sequence<seqan3::dna4>(sequence_length, 0, 1);
    int score = 0;

    for (auto _ : state)
    {
        for (auto && rng : align_pairwise(std::tie(seq1, seq2), edit_distance_cfg))
            score += rng.score();
    }

    state.counters["score"] = score;
    state.counters["cells"] = pairwise_cell_updates(ranges::view::single(std::tie(seq1, seq2)), edit_distance_cfg);
    state.counters["CUPS"] = cell_updates_per_second(state.counters["cells"]);
}

#ifdef SEQAN3_HAS_SEQAN2
void seqan2_edit_distance_dna4(benchmark::State & state)
{
    size_t sequence_length = 500;
    auto seq1 = generate_sequence_seqan2<seqan::Dna>(sequence_length, 0, 0);
    auto seq2 = generate_sequence_seqan2<seqan::Dna>(sequence_length, 0, 1);
    int score = 0;

    for (auto _ : state)
        score += global_edit_distance_seqan2(seq1, seq2);

    state.counters["score"] = score;
    state.counters["cells"] = pairwise_cell_updates(ranges::view::single(std::tie(seq1, seq2)), edit_distance_cfg);
    state.counters["CUPS"] = cell_updates_per_second(state.counters["cells"]);
}

void seqan2_edit_distance_generic_dna4(benchmark::State & state)
{
    size_t sequence_length = 500;
    auto seq1 = generate_sequence_seqan2<seqan::Dna>(sequence_length, 0, 0);
    auto seq2 = generate_sequence_seqan2<seqan::Dna>(sequence_length, 0, 1);
    int score = 0;

    for (auto _ : state)
        score += seqan::globalAlignmentScore(seq1, seq2, seqan::Score<int>{0, -1, -1});

    state.counters["score"] = score;
    state.counters["cells"] = pairwise_cell_updates(ranges::view::single(std::tie(seq1, seq2)), edit_distance_cfg);
    state.counters["CUPS"] = cell_updates_per_second(state.counters["cells"]);
}
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  edit_distance; score; dna4; set
// ============================================================================

void seqan3_edit_distance_dna4_collection(benchmark::State & state)
{
    size_t sequence_length = 500;
    size_t set_size = 100;

    auto vec = generate_sequence_pairs<seqan3::dna4>(sequence_length, set_size);
    int score = 0;

    for (auto _ : state)
    {
        for (auto && [seq1, seq2] : vec)
        {
            detail::edit_distance_unbanded edit_distance{seq1, seq2, edit_distance_cfg};
            edit_distance(0u);
            score += edit_distance.score().value_or(0u);
        }
    }

    state.counters["score"] = score;
    state.counters["cells"] = pairwise_cell_updates(vec, edit_distance_cfg);
    state.counters["CUPS"] = cell_updates_per_second(state.counters["cells"]);
}

void seqan3_edit_distance_dna4_collection_selector(benchmark::State & state)
{
    size_t sequence_length = 500;
    size_t set_size = 100;

    auto vec = generate_sequence_pairs<seqan3::dna4>(sequence_length, set_size);
    int score = 0;

    for (auto _ : state)
    {
        for (auto && rng : align_pairwise(vec, edit_distance_cfg))
            score += rng.score();
    }

    state.counters["score"] = score;
    state.counters["cells"] = pairwise_cell_updates(vec, edit_distance_cfg);
    state.counters["CUPS"] = cell_updates_per_second(state.counters["cells"]);
}

#ifdef SEQAN3_HAS_SEQAN2
void seqan2_edit_distance_dna4_collection(benchmark::State & state)
{
    size_t sequence_length = 500;
    size_t set_size = 100;

    auto [vec1, vec2] = generate_sequence_pairs_seqan2<seqan::Dna>(sequence_length, set_size);
    int score = 0;

    for (auto _ : state)
    {
        for (unsigned i = 0; i < set_size; ++i)
            score += global_edit_distance_seqan2(vec1[i], vec2[i]);
    }

    state.counters["score"] = score;
    state.counters["cells"] = pairwise_cell_updates(ranges::view::zip(vec1, vec2), edit_distance_cfg);
    state.counters["CUPS"] = cell_updates_per_second(state.counters["cells"]);
}

void seqan2_edit_distance_dna4_generic_collection(benchmark::State & state)
{
    size_t sequence_length = 500;
    size_t set_size = 100;

    auto [vec1, vec2] = generate_sequence_pairs_seqan2<seqan::Dna>(sequence_length, set_size);
    int score = 0;

    for (auto _ : state)
    {
        for (int score_: seqan::globalAlignmentScore(vec1, vec2, seqan::Score<int>{0, -1, -1}))
            score += score_;
    }

    state.counters["score"] = score;
    state.counters["cells"] = pairwise_cell_updates(ranges::view::zip(vec1, vec2), edit_distance_cfg);
    state.counters["CUPS"] = cell_updates_per_second(state.counters["cells"]);
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
