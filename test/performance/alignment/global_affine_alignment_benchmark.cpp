// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <iostream>
#include <memory>
#include <random>
#include <ranges>
#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/utility/views/zip.hpp>

#ifdef SEQAN3_HAS_SEQAN2
#    include <seqan/align.h>
#endif

constexpr auto nt_score_scheme = seqan3::nucleotide_scoring_scheme{seqan3::match_score{4}, seqan3::mismatch_score{-5}};
constexpr auto affine_cfg =
    seqan3::align_cfg::method_global{}
    | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10}, seqan3::align_cfg::extension_score{-1}}
    | seqan3::align_cfg::scoring_scheme{nt_score_scheme};

// ============================================================================
//  affine; score; dna4; single
// ============================================================================

void seqan3_affine_dna4(benchmark::State & state)
{
    size_t sequence_length = 500;
    auto seq1 = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 0);
    auto seq2 = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 1);

    for (auto _ : state)
    {
        auto rng = align_pairwise(std::tie(seq1, seq2), affine_cfg | seqan3::align_cfg::output_score{});
        *std::ranges::begin(rng);
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(std::views::single(std::tie(seq1, seq2)), affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(seqan3_affine_dna4);

#ifdef SEQAN3_HAS_SEQAN2

void seqan2_affine_dna4(benchmark::State & state)
{
    size_t sequence_length = 500;
    auto seq1 = seqan3::test::generate_sequence_seqan2<seqan2::Dna>(sequence_length, 0, 0);
    auto seq2 = seqan3::test::generate_sequence_seqan2<seqan2::Dna>(sequence_length, 0, 1);

    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan2::globalAlignmentScore(seq1, seq2, seqan2::Score<int>{4, -5, -1, -11});
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(std::views::single(std::tie(seq1, seq2)), affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(seqan2_affine_dna4);
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  affine; trace; dna4; single
// ============================================================================

void seqan3_affine_dna4_trace(benchmark::State & state)
{
    size_t sequence_length = 500;
    auto seq1 = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 0);
    auto seq2 = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 1);

    for (auto _ : state)
    {
        auto rng = align_pairwise(std::tie(seq1, seq2), affine_cfg | seqan3::align_cfg::output_alignment{});
        *std::ranges::begin(rng);
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(std::views::single(std::tie(seq1, seq2)), affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(seqan3_affine_dna4_trace);

#ifdef SEQAN3_HAS_SEQAN2

void seqan2_affine_dna4_trace(benchmark::State & state)
{
    size_t sequence_length = 500;
    auto seq1 = seqan3::test::generate_sequence_seqan2<seqan2::Dna>(sequence_length, 0, 0);
    auto seq2 = seqan3::test::generate_sequence_seqan2<seqan2::Dna>(sequence_length, 0, 1);

    seqan2::Gaps<decltype(seq1)> gap1{seq1};
    seqan2::Gaps<decltype(seq2)> gap2{seq2};
    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan2::globalAlignment(gap1, gap2, seqan2::Score<int>{4, -5, -1, -11});
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(std::views::single(std::tie(seq1, seq2)), affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(seqan2_affine_dna4_trace);
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  affine; score; dna4; collection
// ============================================================================

void seqan3_affine_dna4_collection(benchmark::State & state)
{
    size_t sequence_length = 100;
    size_t set_size = 100;
    using sequence_t = decltype(seqan3::test::generate_sequence<seqan3::dna4>());

    std::vector<std::pair<sequence_t, sequence_t>> vec;
    for (unsigned i = 0; i < set_size; ++i)
    {
        sequence_t seq1 = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, i);
        sequence_t seq2 = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, i + set_size);
        vec.push_back(std::pair{seq1, seq2});
    }

    for (auto _ : state)
    {
        for (auto && rng : align_pairwise(vec, affine_cfg | seqan3::align_cfg::output_score{}))
            rng.score();
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(vec, affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(seqan3_affine_dna4_collection);

#ifdef SEQAN3_HAS_SEQAN2

void seqan2_affine_dna4_collection(benchmark::State & state)
{
    size_t sequence_length = 100;
    size_t set_size = 100;
    using sequence_t = decltype(seqan3::test::generate_sequence_seqan2<seqan2::Dna>());

    seqan2::StringSet<sequence_t> vec1;
    seqan2::StringSet<sequence_t> vec2;
    for (unsigned i = 0; i < set_size; ++i)
    {
        sequence_t seq1 = seqan3::test::generate_sequence_seqan2<seqan2::Dna>(sequence_length, 0, i);
        sequence_t seq2 = seqan3::test::generate_sequence_seqan2<seqan2::Dna>(sequence_length, 0, i + set_size);
        appendValue(vec1, seq1);
        appendValue(vec2, seq2);
    }

    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan2::globalAlignmentScore(vec1, vec2, seqan2::Score<int>{4, -5, -1, -11});
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(seqan3::views::zip(vec1, vec2), affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(seqan2_affine_dna4_collection);
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  affine; trace; dna4; collection
// ============================================================================

void seqan3_affine_dna4_trace_collection(benchmark::State & state)
{
    size_t sequence_length = 100;
    size_t set_size = 100;
    using sequence_t = decltype(seqan3::test::generate_sequence<seqan3::dna4>());

    std::vector<std::pair<sequence_t, sequence_t>> vec;
    for (unsigned i = 0; i < set_size; ++i)
    {
        sequence_t seq1 = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, i);
        sequence_t seq2 = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, i + set_size);
        vec.push_back(std::pair{seq1, seq2});
    }

    for (auto _ : state)
    {
        for (auto && rng : align_pairwise(vec, affine_cfg | seqan3::align_cfg::output_alignment{}))
            rng.alignment();
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(vec, affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(seqan3_affine_dna4_trace_collection);

#ifdef SEQAN3_HAS_SEQAN2

void seqan2_affine_dna4_trace_collection(benchmark::State & state)
{
    size_t sequence_length = 100;
    size_t set_size = 100;
    using sequence_t = decltype(seqan3::test::generate_sequence_seqan2<seqan2::Dna>());

    seqan2::StringSet<sequence_t> vec1;
    seqan2::StringSet<sequence_t> vec2;
    for (unsigned i = 0; i < set_size; ++i)
    {
        sequence_t seq1 = seqan3::test::generate_sequence_seqan2<seqan2::Dna>(sequence_length, 0, i);
        sequence_t seq2 = seqan3::test::generate_sequence_seqan2<seqan2::Dna>(sequence_length, 0, i + set_size);
        appendValue(vec1, seq1);
        appendValue(vec2, seq2);
    }

    seqan2::StringSet<seqan2::Gaps<sequence_t>> gap1;
    seqan2::StringSet<seqan2::Gaps<sequence_t>> gap2;

    for (unsigned i = 0; i < set_size; ++i)
    {
        appendValue(gap1, seqan2::Gaps<sequence_t>{vec1[i]});
        appendValue(gap2, seqan2::Gaps<sequence_t>{vec2[i]});
    }

    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan2::globalAlignment(gap1, gap2, seqan2::Score<int>{4, -5, -1, -11});
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(seqan3::views::zip(vec1, vec2), affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(seqan2_affine_dna4_trace_collection);
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
