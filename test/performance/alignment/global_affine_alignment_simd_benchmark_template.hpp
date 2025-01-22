// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <ranges>
#include <tuple>
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
#    include <seqan/align_parallel.h>
#endif

// Globally defined constants to ensure same test data.
inline constexpr size_t sequence_length = 150;
#ifndef NDEBUG
inline constexpr size_t set_size = 16;
#else
inline constexpr size_t set_size = 1024;
#endif // NDEBUG

// We don't know if the system supports hyper-threading so we use only half the threads so that the
// simd benchmark is likely to run on physical cores only.
uint32_t get_number_of_threads()
{
    uint32_t thread_count = std::thread::hardware_concurrency();
    return (thread_count == 1) ? thread_count : thread_count >> 1;
}

// ----------------------------------------------------------------------------
//  seqan3 pairwise alignment
// ----------------------------------------------------------------------------

template <typename alphabet_t, typename... align_configs_t>
void seqan3_affine_accelerated(benchmark::State & state, alphabet_t, align_configs_t &&... configs)
{
    size_t sequence_length_variance = state.range(0);
    auto data = seqan3::test::generate_sequence_pairs<alphabet_t>(sequence_length, set_size, sequence_length_variance);

    int64_t total = 0;
    auto accelerate_config = (configs | ...);
    for (auto _ : state)
    {
        for (auto && res : seqan3::align_pairwise(data, accelerate_config))
            total += res.score();
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(data, accelerate_config);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

#ifdef SEQAN3_HAS_SEQAN2

// ----------------------------------------------------------------------------
//  seqan2 pairwise alignment
// ----------------------------------------------------------------------------

template <typename alphabet_t, typename... args_t>
void seqan2_affine_accelerated(benchmark::State & state, alphabet_t, args_t &&... args)
{
    std::tuple captured_args{args...};

    size_t sequence_length_variance = state.range(0);
    auto [vec1, vec2] =
        seqan3::test::generate_sequence_pairs_seqan2<alphabet_t>(sequence_length, set_size, sequence_length_variance);

    auto scoring_scheme = std::get<0>(captured_args);
    auto exec = std::get<1>(captured_args);
    setNumThreads(exec, std::get<2>(captured_args));

    // Possibly enable banded alignment.
    auto seqan3_align_cfg = std::get<3>(captured_args);
    static constexpr bool execute_with_band = seqan3_align_cfg.template exists<seqan3::align_cfg::band_fixed_size>();
    [[maybe_unused]] int lower_diagonal{};
    [[maybe_unused]] int upper_diagonal{};

    if constexpr (execute_with_band)
    {
        using std::get;
        auto band_cfg = get<seqan3::align_cfg::band_fixed_size>(seqan3_align_cfg);
        lower_diagonal = band_cfg.lower_diagonal;
        upper_diagonal = band_cfg.upper_diagonal;
    }

    int64_t total = 0;
    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan2::String<int> res;
        if constexpr (execute_with_band)
            res = seqan2::globalAlignmentScore(exec, vec1, vec2, scoring_scheme, lower_diagonal, upper_diagonal);
        else
            res = seqan2::globalAlignmentScore(exec, vec1, vec2, scoring_scheme);

        total = std::accumulate(seqan2::begin(res), seqan2::end(res), total);
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(seqan3::views::zip(vec1, vec2), seqan3_align_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

#endif // SEQAN3_HAS_SEQAN2
