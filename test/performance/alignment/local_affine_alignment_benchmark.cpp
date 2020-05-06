// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <iostream>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/zip.hpp>

#include <seqan3/test/performance/units.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/std/ranges>

#ifdef SEQAN3_HAS_SEQAN2
    #include <seqan/align.h>
#endif

constexpr auto local_affine_cfg = seqan3::align_cfg::mode{seqan3::local_alignment} |
                                  seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1},
                                                                            seqan3::gap_open_score{-10}}} |
                                  seqan3::align_cfg::scoring{
                                        seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                          seqan3::mismatch_score{-5}}
                                  };

// ============================================================================
//  affine; score; dna4; single
// ============================================================================

void seqan3_affine_dna4(benchmark::State & state)
{
    auto seq1 = seqan3::test::generate_sequence<seqan3::dna4>(500, 0, 0);
    auto seq2 = seqan3::test::generate_sequence<seqan3::dna4>(250, 0, 1);

    for (auto _ : state)
    {
        auto rng = align_pairwise(std::tie(seq1, seq2), local_affine_cfg |
                                                        seqan3::align_cfg::result{seqan3::with_score});
        *std::ranges::begin(rng);
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(std::views::single(std::tie(seq1, seq2)),
                                                                  local_affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(seqan3_affine_dna4);

#ifdef SEQAN3_HAS_SEQAN2

void seqan2_affine_dna4(benchmark::State & state)
{
    auto seq1 = generate_sequence_seqan2<seqan::Dna>(500, 0, 0);
    auto seq2 = generate_sequence_seqan2<seqan::Dna>(250, 0, 1);

    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan::localAlignmentScore(seq1, seq2, seqan::Score<int>{4, -5, -1, -11});
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(std::views::single(std::tie(seq1, seq2)),
                                                                  local_affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(seqan2_affine_dna4);
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  affine; trace; dna4; single
// ============================================================================

void seqan3_affine_dna4_trace(benchmark::State & state)
{
    auto seq1 = seqan3::test::generate_sequence<seqan3::dna4>(500, 0, 0);
    auto seq2 = seqan3::test::generate_sequence<seqan3::dna4>(250, 0, 1);

    for (auto _ : state)
    {
        auto rng = align_pairwise(std::tie(seq1, seq2), local_affine_cfg |
                                                        seqan3::align_cfg::result{seqan3::with_alignment});
        *std::ranges::begin(rng);
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(std::views::single(std::tie(seq1, seq2)),
                                                                  local_affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(seqan3_affine_dna4_trace);

#ifdef SEQAN3_HAS_SEQAN2

void seqan2_affine_dna4_trace(benchmark::State & state)
{
    auto seq1 = generate_sequence_seqan2<seqan::Dna>(500, 0, 0);
    auto seq2 = generate_sequence_seqan2<seqan::Dna>(250, 0, 1);

    seqan::Gaps<decltype(seq1)> gap1{seq1};
    seqan::Gaps<decltype(seq2)> gap2{seq2};
    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan::localAlignment(gap1, gap2, seqan::Score<int>{4, -5, -1, -11});
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(std::views::single(std::tie(seq1, seq2)),
                                                                  local_affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(seqan2_affine_dna4_trace);
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  affine; score; dna4; collection
// ============================================================================

void seqan3_affine_dna4_collection(benchmark::State & state)
{
    using sequence_t = decltype(seqan3::test::generate_sequence<seqan3::dna4>());

    std::vector<std::pair<sequence_t, sequence_t>> vec;
    for (unsigned i = 0; i < 100; ++i)
    {
        sequence_t seq1 = seqan3::test::generate_sequence<seqan3::dna4>(100, 0, i);
        sequence_t seq2 = seqan3::test::generate_sequence<seqan3::dna4>(50, 0, i + 100);
        vec.push_back(std::pair{seq1, seq2});
    }

    for (auto _ : state)
    {
        for (auto && rng : align_pairwise(vec, local_affine_cfg | seqan3::align_cfg::result{seqan3::with_score}))
            rng.score();
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(vec, local_affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(seqan3_affine_dna4_collection);

#ifdef SEQAN3_HAS_SEQAN2

void seqan2_affine_dna4_collection(benchmark::State & state)
{
    using sequence_t = decltype(generate_sequence_seqan2<seqan::Dna>());

    seqan::StringSet<sequence_t> vec1;
    seqan::StringSet<sequence_t> vec2;
    for (unsigned i = 0; i < 100; ++i)
    {
        sequence_t seq1 = generate_sequence_seqan2<seqan::Dna>(100, 0, i);
        sequence_t seq2 = generate_sequence_seqan2<seqan::Dna>(50, 0, i + 100);
        appendValue(vec1, seq1);
        appendValue(vec2, seq2);
    }

    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan::localAlignmentScore(vec1, vec2, seqan::Score<int>{4, -5, -1, -11});
    }
}

BENCHMARK(seqan2_affine_dna4_collection);
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  affine; trace; dna4; collection
// ============================================================================

void seqan3_affine_dna4_trace_collection(benchmark::State & state)
{
    using sequence_t = decltype(seqan3::test::generate_sequence<seqan3::dna4>());

    std::vector<std::pair<sequence_t, sequence_t>> vec;
    for (unsigned i = 0; i < 100; ++i)
    {
        sequence_t seq1 = seqan3::test::generate_sequence<seqan3::dna4>(100, 0, i);
        sequence_t seq2 = seqan3::test::generate_sequence<seqan3::dna4>(50, 0, i + 100);
        vec.push_back(std::pair{seq1, seq2});
    }

    for (auto _ : state)
    {
        for (auto && rng : align_pairwise(vec, local_affine_cfg | seqan3::align_cfg::result{seqan3::with_alignment}))
            rng.score();
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(vec, local_affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(seqan3_affine_dna4_trace_collection);

#ifdef SEQAN3_HAS_SEQAN2

void seqan2_affine_dna4_trace_collection(benchmark::State & state)
{
    using sequence_t = decltype(generate_sequence_seqan2<seqan::Dna>());

    seqan::StringSet<sequence_t> vec1;
    seqan::StringSet<sequence_t> vec2;
    for (unsigned i = 0; i < 100; ++i)
    {
        sequence_t seq1 = generate_sequence_seqan2<seqan::Dna>(100, 0, i);
        sequence_t seq2 = generate_sequence_seqan2<seqan::Dna>(50, 0, i + 100);
        appendValue(vec1, seq1);
        appendValue(vec2, seq2);
    }


    seqan::StringSet<seqan::Gaps<sequence_t>> gap1;
    seqan::StringSet<seqan::Gaps<sequence_t>> gap2;

    for (unsigned i = 0; i < 100; ++i)
    {
        appendValue(gap1, seqan::Gaps<sequence_t>{vec1[i]});
        appendValue(gap2, seqan::Gaps<sequence_t>{vec2[i]});
    }

    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan::localAlignment(gap1, gap2, seqan::Score<int>{4, -5, -1, -11}, seqan::Gotoh());
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(seqan3::views::zip(vec1, vec2), local_affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
}

BENCHMARK(seqan2_affine_dna4_trace_collection);
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
