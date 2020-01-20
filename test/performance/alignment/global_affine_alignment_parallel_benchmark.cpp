// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/std/ranges>

#include <seqan3/core/debug_stream.hpp>

#ifdef SEQAN3_HAS_SEQAN2
    #include <seqan/align.h>
    #include <seqan/align_parallel.h>
#endif

using namespace seqan3;
using namespace seqan3::test;

constexpr auto affine_cfg = align_cfg::mode{global_alignment} |
                            align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
                            align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}};

// Aliases to beautify the benchmark output
using score = detail::with_score_type;
using trace = detail::with_alignment_type;

// Globally defined constants to ensure same test data.
inline constexpr size_t sequence_length = 100;
inline constexpr size_t set_size        = 500;
inline constexpr size_t variance        = 10;

template <typename alphabet_t>
auto generate_data_seqan3()
{
    using sequence_t = decltype(generate_sequence<alphabet_t>());

    std::vector<sequence_t> vec1;
    std::vector<sequence_t> vec2;
    for (unsigned i = 0; i < set_size; ++i)
    {
        vec1.push_back(generate_sequence<alphabet_t>(sequence_length, variance, i));
        vec2.push_back(generate_sequence<alphabet_t>(sequence_length, variance, i + set_size));
    }
    return std::pair{vec1, vec2};
}

#ifdef SEQAN3_HAS_SEQAN2
template <typename alphabet_t>
auto generate_data_seqan2()
{
    using sequence_t = decltype(generate_sequence_seqan2<alphabet_t>());

    seqan::StringSet<sequence_t> vec1;
    seqan::StringSet<sequence_t> vec2;
    for (unsigned i = 0; i < set_size; ++i)
    {
        appendValue(vec1, generate_sequence_seqan2<alphabet_t>(sequence_length, variance, i));
        appendValue(vec2, generate_sequence_seqan2<alphabet_t>(sequence_length, variance, i + set_size));
    }
    return std::pair{vec1, vec2};
}
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  affine; score; dna4; collection
// ============================================================================

template <typename result_t>
void seqan3_affine_dna4_parallel(benchmark::State & state)
{
    auto [vec1, vec2] = generate_data_seqan3<seqan3::dna4>();

    auto data = views::zip(vec1, vec2) | views::to<std::vector>;

    int64_t total = 0;
    for (auto _ : state)
    {
        for (auto && res : align_pairwise(data, affine_cfg |
                                          align_cfg::result{result_t{}} |
                                          align_cfg::parallel{std::thread::hardware_concurrency()}))
        {
            total += res.score();
        }
    }

    state.counters["cells"] = pairwise_cell_updates(views::zip(vec1, vec2), affine_cfg);
    state.counters["CUPS"] = cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

BENCHMARK_TEMPLATE(seqan3_affine_dna4_parallel, score)->UseRealTime();
BENCHMARK_TEMPLATE(seqan3_affine_dna4_parallel, trace)->UseRealTime();

#if defined(_OPENMP)
template <typename result_t>
void seqan3_affine_dna4_omp_for(benchmark::State & state)
{
    auto [vec1, vec2] = generate_data_seqan3<seqan3::dna4>();
    auto zip = views::zip(vec1, vec2);
    int64_t total = 0;
    for (auto _ : state)
    {
        #pragma omp parallel for num_threads(std::thread::hardware_concurrency()) schedule(guided)
        for (size_t i = 0; i < zip.size(); ++i)
        {
            auto rng = align_pairwise(zip[i], affine_cfg | align_cfg::result{result_t{}});
            auto res = *rng.begin();
            total += res.score();
        }
    }

    state.counters["cells"] = pairwise_cell_updates(views::zip(vec1, vec2), affine_cfg);
    state.counters["CUPS"] = cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

BENCHMARK_TEMPLATE(seqan3_affine_dna4_omp_for, score)->UseRealTime();
BENCHMARK_TEMPLATE(seqan3_affine_dna4_omp_for, trace)->UseRealTime();
#endif // defined(_OPENMP)

#ifdef SEQAN3_HAS_SEQAN2

template <typename result_t>
void seqan2_affine_dna4_parallel(benchmark::State & state)
{
    using namespace seqan;

    auto [vec1, vec2] = generate_data_seqan2<seqan::Dna>();
    ExecutionPolicy<Parallel, Serial> exec{};
    setNumThreads(exec, std::thread::hardware_concurrency());

    int64_t total = 0;
    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        auto res = seqan::globalAlignmentScore(exec, vec1, vec2, seqan::Score<int>{4, -5, -1, -11});
        total = std::accumulate(seqan::begin(res), seqan::end(res), total);
    }

    state.counters["cells"] = pairwise_cell_updates(views::zip(vec1, vec2), affine_cfg);
    state.counters["CUPS"] = cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

// Note SeqAn2 has no parallel interface yet for computing the traceback as well.
BENCHMARK_TEMPLATE(seqan2_affine_dna4_parallel, score)->UseRealTime();

#if defined(_OPENMP)
template <typename result_t>
void seqan2_affine_dna4_omp_for(benchmark::State & state)
{
    constexpr bool score_only = std::is_same_v<result_t, seqan3::detail::with_score_type>;

    auto [vec1, vec2] = generate_data_seqan2<seqan::Dna>();

    using gapped_t = std::conditional_t<score_only,
                                        std::remove_reference_t<decltype(vec1[0])>,
                                        seqan::Gaps<std::remove_reference_t<decltype(vec1[0])>>>;

    seqan::StringSet<gapped_t> gap1;
    seqan::StringSet<gapped_t> gap2;

    for (size_t i = 0; i < seqan::length(vec1); ++i)
    {
        appendValue(gap1, gapped_t{vec1[i]});
        appendValue(gap2, gapped_t{vec2[i]});
    }

    int64_t total = 0;
    for (auto _ : state)
    {
        #pragma omp parallel for num_threads(std::thread::hardware_concurrency()) schedule(guided)
        for (size_t i = 0; i < seqan::length(vec1); ++i)
        {
            if constexpr (score_only)
                total += seqan::globalAlignmentScore(vec1[i], vec2[i], seqan::Score<int>{4, -5, -1, -11});
            else
                total += seqan::globalAlignment(gap1[i], gap2[i], seqan::Score<int>{4, -5, -1, -11});
        }
    }

    state.counters["cells"] = pairwise_cell_updates(views::zip(vec1, vec2), affine_cfg);
    state.counters["CUPS"] = cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

BENCHMARK_TEMPLATE(seqan2_affine_dna4_omp_for, score)->UseRealTime();
BENCHMARK_TEMPLATE(seqan2_affine_dna4_omp_for, trace)->UseRealTime();

#endif // defined(_OPENMP)
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
