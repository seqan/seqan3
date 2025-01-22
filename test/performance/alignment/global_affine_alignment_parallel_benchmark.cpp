// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <random>
#include <ranges>
#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/utility/range/to.hpp>
#include <seqan3/utility/views/zip.hpp>

#ifdef SEQAN3_HAS_SEQAN2
#    include <seqan/align.h>
#    include <seqan/align_parallel.h>
#endif

constexpr auto nt_score_scheme = seqan3::nucleotide_scoring_scheme{seqan3::match_score{4}, seqan3::mismatch_score{-5}};
constexpr auto affine_cfg =
    seqan3::align_cfg::method_global{}
    | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10}, seqan3::align_cfg::extension_score{-1}}
    | seqan3::align_cfg::scoring_scheme{nt_score_scheme};

// Aliases to beautify the benchmark output
using score = seqan3::align_cfg::output_score;
using trace = decltype(seqan3::align_cfg::output_score{} | seqan3::align_cfg::output_alignment{});

// Globally defined constants to ensure same test data.
inline constexpr size_t sequence_length = 100;
inline constexpr size_t set_size = 500;
inline constexpr size_t variance = 10;

template <typename alphabet_t>
auto generate_data_seqan3()
{
    using sequence_t = decltype(seqan3::test::generate_sequence<alphabet_t>());

    std::vector<sequence_t> vec1;
    std::vector<sequence_t> vec2;
    for (unsigned i = 0; i < set_size; ++i)
    {
        vec1.push_back(seqan3::test::generate_sequence<alphabet_t>(sequence_length, variance, i));
        vec2.push_back(seqan3::test::generate_sequence<alphabet_t>(sequence_length, variance, i + set_size));
    }
    return std::pair{vec1, vec2};
}

#ifdef SEQAN3_HAS_SEQAN2
template <typename alphabet_t>
auto generate_data_seqan2()
{
    using sequence_t = decltype(seqan3::test::generate_sequence_seqan2<alphabet_t>());

    seqan2::StringSet<sequence_t> vec1;
    seqan2::StringSet<sequence_t> vec2;
    for (unsigned i = 0; i < set_size; ++i)
    {
        appendValue(vec1, seqan3::test::generate_sequence_seqan2<alphabet_t>(sequence_length, variance, i));
        appendValue(vec2, seqan3::test::generate_sequence_seqan2<alphabet_t>(sequence_length, variance, i + set_size));
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

    auto data = seqan3::views::zip(vec1, vec2) | seqan3::ranges::to<std::vector>();

    int64_t total = 0;
    for (auto _ : state)
    {
        for (auto && res :
             align_pairwise(data,
                            affine_cfg | result_t{} | seqan3::align_cfg::parallel{std::thread::hardware_concurrency()}))
        {
            total += res.score();
        }
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(seqan3::views::zip(vec1, vec2), affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

BENCHMARK_TEMPLATE(seqan3_affine_dna4_parallel, score)->UseRealTime();
BENCHMARK_TEMPLATE(seqan3_affine_dna4_parallel, trace)->UseRealTime();

#if defined(_OPENMP)
template <typename result_t>
void seqan3_affine_dna4_omp_for(benchmark::State & state)
{
    auto [vec1, vec2] = generate_data_seqan3<seqan3::dna4>();
    auto zip = seqan3::views::zip(vec1, vec2);
    int64_t total = 0;
    for (auto _ : state)
    {
#    pragma omp parallel for num_threads(std::thread::hardware_concurrency()) schedule(guided) reduction(+:total)
        for (size_t i = 0; i < zip.size(); ++i)
        {
            auto rng = align_pairwise(zip[i], affine_cfg | result_t{});
            auto res = *rng.begin();
            total += res.score();
        }
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(seqan3::views::zip(vec1, vec2), affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

BENCHMARK_TEMPLATE(seqan3_affine_dna4_omp_for, score)->UseRealTime();
BENCHMARK_TEMPLATE(seqan3_affine_dna4_omp_for, trace)->UseRealTime();
#endif // defined(_OPENMP)

#ifdef SEQAN3_HAS_SEQAN2

template <typename result_t>
void seqan2_affine_dna4_parallel(benchmark::State & state)
{
    auto [vec1, vec2] = generate_data_seqan2<seqan2::Dna>();
    seqan2::ExecutionPolicy<seqan2::Parallel, seqan2::Serial> exec{};
    setNumThreads(exec, std::thread::hardware_concurrency());

    int64_t total = 0;
    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        auto res = seqan2::globalAlignmentScore(exec, vec1, vec2, seqan2::Score<int>{4, -5, -1, -11});
        total = std::accumulate(seqan2::begin(res), seqan2::end(res), total);
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(seqan3::views::zip(vec1, vec2), affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

// Note SeqAn2 has no parallel interface yet for computing the traceback as well.
BENCHMARK_TEMPLATE(seqan2_affine_dna4_parallel, score)->UseRealTime();
#endif // SEQAN3_HAS_SEQAN2

// Crashes with libc++ or IntelLLVM
#if defined(SEQAN3_HAS_SEQAN2) && defined(_OPENMP) && !defined(_LIBCPP_VERSION) && !defined(__INTEL_LLVM_COMPILER)
template <typename result_t>
void seqan2_affine_dna4_omp_for(benchmark::State & state)
{
    constexpr bool score_only = std::is_same_v<result_t, seqan3::align_cfg::output_score>;

    auto [vec1, vec2] = generate_data_seqan2<seqan2::Dna>();

    using gapped_t = std::conditional_t<score_only,
                                        std::remove_reference_t<decltype(vec1[0])>,
                                        seqan2::Gaps<std::remove_reference_t<decltype(vec1[0])>>>;

    seqan2::StringSet<gapped_t> gap1;
    seqan2::StringSet<gapped_t> gap2;

    for (size_t i = 0; i < seqan2::length(vec1); ++i)
    {
        appendValue(gap1, gapped_t{vec1[i]});
        appendValue(gap2, gapped_t{vec2[i]});
    }

    int64_t total = 0;
    for (auto _ : state)
    {
#    pragma omp parallel for num_threads(std::thread::hardware_concurrency()) schedule(guided) reduction(+:total)
        for (size_t i = 0; i < seqan2::length(vec1); ++i)
        {
            if constexpr (score_only)
                total += seqan2::globalAlignmentScore(vec1[i], vec2[i], seqan2::Score<int>{4, -5, -1, -11});
            else
                total += seqan2::globalAlignment(gap1[i], gap2[i], seqan2::Score<int>{4, -5, -1, -11});
        }
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(seqan3::views::zip(vec1, vec2), affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

BENCHMARK_TEMPLATE(seqan2_affine_dna4_omp_for, score)->UseRealTime();
BENCHMARK_TEMPLATE(seqan2_affine_dna4_omp_for, trace)->UseRealTime();
#endif // defined(SEQAN3_HAS_SEQAN2) && defined(_OPENMP) && !defined(_LIBCPP_VERSION) && !defined(__INTEL_LLVM_COMPILER)

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
