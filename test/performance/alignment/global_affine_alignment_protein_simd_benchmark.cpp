// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include "global_affine_alignment_simd_benchmark_template.hpp"

// Range to test for sequence length variance
inline constexpr size_t deviation_begin = 0;
inline constexpr size_t deviation_end = 64;
inline constexpr size_t deviation_step = 8;

// ----------------------------------------------------------------------------
// SeqAn3
// ----------------------------------------------------------------------------

constexpr auto aa_score_scheme = seqan3::aminoacid_scoring_scheme{seqan3::aminoacid_similarity_matrix::blosum62};
constexpr auto affine_cfg = seqan3::align_cfg::method_global{} |
                            seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10},
                                                               seqan3::align_cfg::extension_score{-1}} |
                            seqan3::align_cfg::scoring_scheme{aa_score_scheme};

BENCHMARK_CAPTURE(seqan3_affine_accelerated,
                  simd_with_score,
                  seqan3::aa27{},
                  affine_cfg,
                  seqan3::align_cfg::output_score{},
                  seqan3::align_cfg::score_type<int16_t>{},
                  seqan3::align_cfg::vectorised{})
                        ->UseRealTime()
                        ->DenseRange(deviation_begin, deviation_end, deviation_step);

BENCHMARK_CAPTURE(seqan3_affine_accelerated,
                  simd_parallel_with_score,
                  seqan3::aa27{},
                  affine_cfg,
                  seqan3::align_cfg::output_score{},
                  seqan3::align_cfg::vectorised{},
                  seqan3::align_cfg::score_type<int16_t>{},
                  seqan3::align_cfg::parallel{get_number_of_threads()})
                        ->UseRealTime()
                        ->DenseRange(deviation_begin, deviation_end, deviation_step);

#ifdef SEQAN3_HAS_SEQAN2

// ----------------------------------------------------------------------------
// SeqAn2
// ----------------------------------------------------------------------------

using scoring_scheme_t = seqan::Score<int16_t, seqan::ScoreMatrix<seqan::AminoAcid, seqan::ScoreSpecBlosum62>>;

// Note SeqAn2 has no parallel interface yet for computing the traceback as well.
BENCHMARK_CAPTURE(seqan2_affine_accelerated,
                  simd_with_score,
                  seqan::AminoAcid{},
                  scoring_scheme_t{-1, -11},
                  seqan::ExecutionPolicy<seqan::Serial, seqan::Vectorial>{},
                  1,
                  affine_cfg)
                        ->UseRealTime()
                        ->DenseRange(deviation_begin, deviation_end, deviation_step);

BENCHMARK_CAPTURE(seqan2_affine_accelerated,
                  simd_parallel_with_score,
                  seqan::AminoAcid{},
                  scoring_scheme_t{-1, -11},
                  seqan::ExecutionPolicy<seqan::Parallel, seqan::Vectorial>{},
                  get_number_of_threads(),
                  affine_cfg)
                        ->UseRealTime()
                        ->DenseRange(deviation_begin, deviation_end, deviation_step);
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
