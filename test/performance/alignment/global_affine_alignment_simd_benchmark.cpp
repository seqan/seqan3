// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <tuple>
#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/std/ranges>

#ifdef SEQAN3_HAS_SEQAN2
    #include <seqan/align.h>
    #include <seqan/align_parallel.h>
#endif

constexpr auto affine_cfg = seqan3::align_cfg::mode{seqan3::global_alignment} |
                            seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1},
                                                                      seqan3::gap_open_score{-10}}} |
                            seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                         seqan3::mismatch_score{-5}}};

// Globally defined constants to ensure same test data.
inline constexpr size_t sequence_length = 150;
inline constexpr size_t set_size        = 1024;

// Range to test for sequence length variance
inline constexpr size_t deviation_begin = 0;
inline constexpr size_t deviation_end = 64;
inline constexpr size_t deviation_step = 8;

//We don't know if the system supports hyper-threading so we use only half the threads so that the
// simd benchmark is likely to run on physical cores only.
uint32_t get_number_of_threads()
{
    uint32_t thread_count = std::thread::hardware_concurrency();
    return (thread_count == 1) ? thread_count : thread_count >> 1;
}

// ============================================================================
//  affine; score; dna4; collection
// ============================================================================

template <typename ...align_configs_t>
void seqan3_affine_dna4_accelerated(benchmark::State & state, align_configs_t && ...configs)
{
    size_t sequence_length_variance = state.range(0);
    auto data = seqan3::test::generate_sequence_pairs<seqan3::dna4>(sequence_length,
                                                                    set_size,
                                                                    sequence_length_variance);

    int64_t total = 0;
    auto accelerate_config = (affine_cfg | ... | configs);
    for (auto _ : state)
    {
        for (auto && res : seqan3::align_pairwise(data, accelerate_config))
            total += res.score();
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(data, affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

BENCHMARK_CAPTURE(seqan3_affine_dna4_accelerated,
                  simd_with_score,
                  seqan3::align_cfg::result{seqan3::with_score, seqan3::using_score_type<int16_t>},
                  seqan3::align_cfg::vectorise)
                        ->UseRealTime()
                        ->DenseRange(deviation_begin, deviation_end, deviation_step);

BENCHMARK_CAPTURE(seqan3_affine_dna4_accelerated,
                  simd_with_end_position,
                  seqan3::align_cfg::result{seqan3::with_back_coordinate, seqan3::using_score_type<int16_t>},
                  seqan3::align_cfg::vectorise)
                        ->UseRealTime()
                        ->DenseRange(deviation_begin, deviation_end, deviation_step);

BENCHMARK_CAPTURE(seqan3_affine_dna4_accelerated,
                  simd_parallel_with_score,
                  seqan3::align_cfg::result{seqan3::with_score, seqan3::using_score_type<int16_t>},
                  seqan3::align_cfg::vectorise,
                  seqan3::align_cfg::parallel{get_number_of_threads()})
                        ->UseRealTime()
                        ->DenseRange(deviation_begin, deviation_end, deviation_step);

BENCHMARK_CAPTURE(seqan3_affine_dna4_accelerated,
                  simd_parallel_with_end_position,
                  seqan3::align_cfg::result{seqan3::with_back_coordinate, seqan3::using_score_type<int16_t>},
                  seqan3::align_cfg::vectorise,
                  seqan3::align_cfg::parallel{get_number_of_threads()})
                        ->UseRealTime()
                        ->DenseRange(deviation_begin, deviation_end, deviation_step);

#ifdef SEQAN3_HAS_SEQAN2

template <typename ...args_t>
void seqan2_affine_dna4_accelerated(benchmark::State & state, args_t && ...args)
{
    std::tuple captured_args{args...};

    size_t sequence_length_variance = state.range(0);
    auto [vec1, vec2] = seqan3::test::generate_sequence_pairs_seqan2<seqan::Dna>(sequence_length,
                                                                                 set_size,
                                                                                 sequence_length_variance);
    auto exec = std::get<0>(captured_args);
    setNumThreads(exec, std::get<1>(captured_args));

    int64_t total = 0;
    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        auto res = seqan::globalAlignmentScore(exec, vec1, vec2, seqan::Score<int>{4, -5, -1, -11});
        total = std::accumulate(seqan::begin(res), seqan::end(res), total);
    }

    state.counters["cells"] = seqan3::test::pairwise_cell_updates(seqan3::views::zip(vec1, vec2), affine_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

// Note SeqAn2 has no parallel interface yet for computing the traceback as well.
BENCHMARK_CAPTURE(seqan2_affine_dna4_accelerated,
                  simd_with_score,
                  seqan::ExecutionPolicy<seqan::Serial, seqan::Vectorial>{},
                  1)
                        ->UseRealTime()
                        ->DenseRange(deviation_begin, deviation_end, deviation_step);

BENCHMARK_CAPTURE(seqan2_affine_dna4_accelerated,
                  simd_parallel_with_score,
                  seqan::ExecutionPolicy<seqan::Parallel, seqan::Vectorial>{},
                  get_number_of_threads())
                        ->UseRealTime()
                        ->DenseRange(deviation_begin, deviation_end, deviation_step);
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
