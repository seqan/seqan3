// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "alignment_fixture.hpp"

using seqan3::operator""_dna4;

namespace seqan3::test::alignment::fixture::global::affine::banded
{
// clang-format off

inline constexpr auto align_config = seqan3::align_cfg::method_global{} |
                                     seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10},
                                                                        seqan3::align_cfg::extension_score{-1}};
inline constexpr auto nt_score_scheme = seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                          seqan3::mismatch_score{-5}};

static auto dna4_01 = []()
{   //    AACCGGTTAACCGGTT
    //   01234567890123456|
    //  0        x        |
    // A1         x       |
    // C2          x      |
    // G3x          x     |
    // T4 x          x    |
    // A5  x          x   |
    // C6   x          x  |
    // G7    x          x |
    // T8     x          x|
    // A9      x          |
    return alignment_fixture
    {
        "AACCGGTTAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        align_config | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-3},
                                                          seqan3::align_cfg::upper_diagonal{8}} |
                       seqan3::align_cfg::scoring_scheme{nt_score_scheme},
        -18,
        "A---ACCGGTTAACCGGTT",
        "ACGTAC----------GTA",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 16u,
        /*.sequence2_end_position = */ 9u,
        std::vector<std::optional<int32_t>>
        {
        //   e  ,A  ,A  ,C  ,C  ,G  ,G  ,T  ,T  ,A  ,A  ,C  ,C  ,G  ,G  ,T  ,T  ,
        /*e*/0  ,-11,-12,-13,-14,-15,-16,-17,-18,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/-11,4  ,-7 ,-8 ,-9 ,-10,-11,-12,-13,-14,INF,INF,INF,INF,INF,INF,INF,
        /*C*/-12,-7 ,-1 ,-3 ,-4 ,-14,-15,-16,-17,-18,-19,INF,INF,INF,INF,INF,INF,
        /*G*/-13,-8 ,-12,-6 ,-8 ,0  ,-10,-12,-13,-14,-15,-16,INF,INF,INF,INF,INF,
        /*T*/INF,-9 ,-13,-15,-11,-11,-5 ,-6 ,-8 ,-18,-19,-20,-21,INF,INF,INF,INF,
        /*A*/INF,INF,-5 ,-16,-17,-12,-16,-10,-11,-4 ,-14,-16,-17,-18,INF,INF,INF,
        /*C*/INF,INF,INF,-1 ,-12,-13,-14,-15,-15,-15,-9 ,-10,-12,-21,-22,INF,INF,
        /*G*/INF,INF,INF,INF,-6 ,-8 ,-9 ,-19,-20,-16,-20,-14,-15,-8 ,-17,-20,INF,
        /*T*/INF,INF,INF,INF,INF,-11,-13,-5 ,-15,-17,-18,-19,-19,-19,-13,-13,-16,
        /*A*/INF,INF,INF,INF,INF,INF,-16,-16,-10,-11,-13,-23,-24,-20,-24,-18,-18
        },
        std::vector<std::optional<seqan3::detail::trace_directions>>
        {
        //   e  ,A  ,A  ,C  ,C  ,G  ,G  ,T  ,T  ,A  ,A  ,C  ,C  ,G  ,G  ,T  ,T  ,
        /*e*/N  ,L  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/U  ,DUL,DUL,l  ,l  ,l  ,l  ,l  ,l  ,D  ,INF,INF,INF,INF,INF,INF,INF,
        /*C*/u  ,UL ,DUL,DUL,DUl,DUl,DUl,DUl,DUl,DUl,D  ,INF,INF,INF,INF,INF,INF,
        /*G*/u  ,uL ,DUL,DUl,DUL,Dul,DuL,l  ,l  ,l  ,l  ,l  ,INF,INF,INF,INF,INF,
        /*T*/INF,u  ,DuL,ul ,Dul,UL ,DUL,DUL,DUl,DUl,DUl,DUl,D  ,INF,INF,INF,INF,
        /*A*/INF,INF,Du ,uL ,ul ,ul ,DUl,DUl,DUl,Dul,DuL,l  ,l  ,l  ,INF,INF,INF,
        /*C*/INF,INF,INF,Du ,DuL,ul ,l  ,l  ,Dul,Ul ,DUl,DUl,DUl,l  ,l  ,INF,INF,
        /*G*/INF,INF,INF,INF,Du ,DuL,Dul,Dul,Dul,ul ,DUl,DUl,DUl,Dul,DUL,l  ,INF,
        /*T*/INF,INF,INF,INF,INF,Du ,DuL,Dul,DuL,ul ,l  ,l  ,Dul,Ul ,DUl,DUl,D  ,
        /*A*/INF,INF,INF,INF,INF,INF,Du ,UL ,DuL,DuL,Dul,Dul,Dul,ul ,DUl,DUl,DUl
        }
    };
}();

static auto dna4_same_sequence_upper_diagonal_0 = []()
{
    //   0123456789|
    //  0x         |
    //  1 x        |
    //  2  x       |
    //  3x  x      |
    //  4 x  x     |
    //  5  x  x    |
    //  6   x  x   |
    //  7    x  x  |
    //  8     x  x |
    //  9      x  x|
    return alignment_fixture
    {
        "ACGTACGTA"_dna4,
        "ACGTACGTA"_dna4,
        align_config | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-3},
                                                          seqan3::align_cfg::upper_diagonal{0}} |
                       seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                    seqan3::mismatch_score{-5}}},
        36,
        "ACGTACGTA",
        "ACGTACGTA",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 9u,
        /*.sequence2_end_position = */ 9u,
        std::vector<std::optional<int32_t>>
        {
        //   e  ,A  ,C  ,G  ,T  ,A  ,C  ,G  ,T  ,A  ,
        /*e*/0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/-11,4  ,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/-12,-7 ,8  ,INF,INF,INF,INF,INF,INF,INF,
        /*G*/-13,-8 ,-3 ,12 ,INF,INF,INF,INF,INF,INF,
        /*T*/INF,-9 ,-4 ,1  ,16 ,INF,INF,INF,INF,INF,
        /*A*/INF,INF,-5 ,0  ,5  ,20 ,INF,INF,INF,INF,
        /*C*/INF,INF,INF,-1 ,4  ,9  ,24 ,INF,INF,INF,
        /*G*/INF,INF,INF,INF,3  ,8  ,13 ,28 ,INF,INF,
        /*T*/INF,INF,INF,INF,INF,7  ,12 ,17 ,32 ,INF,
        /*A*/INF,INF,INF,INF,INF,INF,11 ,16 ,21 ,36
        },
        std::vector<std::optional<seqan3::detail::trace_directions>>
        {
        //   e  ,A  ,C  ,G  ,T  ,A  ,C  ,G  ,T  ,A  ,
        /*e*/N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/U  ,D  ,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/u  ,UL ,D  ,INF,INF,INF,INF,INF,INF,INF,
        /*G*/u  ,uL ,UL ,D  ,INF,INF,INF,INF,INF,INF,
        /*T*/INF,u  ,uL ,UL ,D  ,INF,INF,INF,INF,INF,
        /*A*/INF,INF,u  ,uL ,UL ,D  ,INF,INF,INF,INF,
        /*C*/INF,INF,INF,u  ,uL ,UL ,D  ,INF,INF,INF,
        /*G*/INF,INF,INF,INF,u  ,uL ,UL ,D  ,INF,INF,
        /*T*/INF,INF,INF,INF,INF,u  ,uL ,UL ,D  ,INF,
        /*A*/INF,INF,INF,INF,INF,INF,u  ,uL ,UL ,D
        }
    };
}();

static auto dna4_same_sequence_lower_diagonal_0 = []()
{
    //   0123456789|
    //  0x       x |
    //  1 x       x|
    //  2  x       |
    //  3   x      |
    //  4    x     |
    //  5     x    |
    //  6      x   |
    //  7       x  |
    //  8        x |
    //  9         x|
    return alignment_fixture
    {
        "ACGTACGTA"_dna4,
        "ACGTACGTA"_dna4,
        align_config | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{0},
                                                          seqan3::align_cfg::upper_diagonal{8}} |
                       seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                    seqan3::mismatch_score{-5}}},
        36,
        "ACGTACGTA",
        "ACGTACGTA",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 9u,
        /*.sequence2_end_position = */ 9u,
        std::vector<std::optional<int32_t>>
        {
        //   e  ,A  ,C  ,G  ,T  ,A  ,C  ,G  ,T  ,A  ,
        /*e*/0  ,-11,-12,-13,-14,-15,-16,-17,-18,INF,
        /*A*/INF,4  ,-7 ,-8 ,-9 ,-10,-11,-12,-13,-14,
        /*C*/INF,INF,8  ,-3 ,-4 ,-5 ,-6 ,-7 ,-8 ,-9 ,
        /*G*/INF,INF,INF,12 ,1  ,0  ,-1 ,-2 ,-3 ,-4 ,
        /*T*/INF,INF,INF,INF,16 ,5  ,4  ,3  ,2  ,1  ,
        /*A*/INF,INF,INF,INF,INF,20 ,9  ,8  ,7  ,6  ,
        /*C*/INF,INF,INF,INF,INF,INF,24 ,13 ,12 ,11 ,
        /*G*/INF,INF,INF,INF,INF,INF,INF,28 ,17 ,16 ,
        /*T*/INF,INF,INF,INF,INF,INF,INF,INF,32 ,21 ,
        /*A*/INF,INF,INF,INF,INF,INF,INF,INF,INF,36
        },
        std::vector<std::optional<seqan3::detail::trace_directions>>
        {
        //   e  ,A  ,C  ,G  ,T  ,A  ,C  ,G  ,T  ,A  ,
        /*e*/N  ,L  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,INF,
        /*A*/INF,DU ,L  ,l  ,l  ,DUl,l  ,l  ,l  ,D  ,
        /*C*/INF,INF,DU ,L  ,l  ,l  ,DUl,l  ,l  ,l  ,
        /*G*/INF,INF,INF,DU ,L  ,l  ,l  ,DUl,l  ,l  ,
        /*T*/INF,INF,INF,INF,DU ,L  ,l  ,l  ,DUl,l  ,
        /*A*/INF,INF,INF,INF,INF,DU ,L  ,l  ,l  ,DUl,
        /*C*/INF,INF,INF,INF,INF,INF,DU ,L  ,l  ,l  ,
        /*G*/INF,INF,INF,INF,INF,INF,INF,DU ,L  ,l  ,
        /*T*/INF,INF,INF,INF,INF,INF,INF,INF,DU ,L  ,
        /*A*/INF,INF,INF,INF,INF,INF,INF,INF,INF,DU
        }
    };
}();

static auto dna4_small_band = []()
{
    //   01234567890|
    // 0   x        |
    // 1    x       |
    // 2x    x      |
    // 3 x    x     |
    // 4  x    x    |
    // 5   x    x   |
    // 6    x    x  |
    // 7     x    x |
    // 8      x    x|
    // 9       x    |
    return alignment_fixture
    {
        "ATCGACGATA"_dna4,
        "ACGACTAGC"_dna4,
        align_config | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-2},
                                                          seqan3::align_cfg::upper_diagonal{2}} |
                       seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                    seqan3::mismatch_score{-5}}},
        -2,
        "ATCGACGATA",
        "A-CGACTAGC",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 10u,
        /*.sequence2_end_position = */ 9u,
        std::vector<std::optional<int32_t>>
        {
        //   e  ,A  ,T  ,C  ,G  ,A  ,C  ,G  ,A  ,T  ,A  ,
        /*e*/0  ,-11,-12,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/-11,4  ,-7 ,-8 ,INF,INF,INF,INF,INF,INF,INF,
        /*C*/-12,-7 ,-1 ,-3 ,-13,INF,INF,INF,INF,INF,INF,
        /*G*/INF,-8 ,-12,-6 ,1  ,-10,INF,INF,INF,INF,INF,
        /*A*/INF,INF,-13,-15,-10,5  ,-6 ,INF,INF,INF,INF,
        /*C*/INF,INF,INF,-9 ,-11,-6 ,9  ,-2 ,INF,INF,INF,
        /*T*/INF,INF,INF,INF,-12,-7 ,-2 ,4  ,-7 ,INF,INF,
        /*A*/INF,INF,INF,INF,INF,-8 ,-3 ,-7 ,8  ,-3 ,INF,
        /*G*/INF,INF,INF,INF,INF,INF,-4 ,1  ,-3 ,3  ,-8 ,
        /*C*/INF,INF,INF,INF,INF,INF,INF,-9 ,-4 ,-8 ,-2
        },
        std::vector<std::optional<seqan3::detail::trace_directions>>
        {
        //   e  ,A  ,T  ,C  ,G  ,A  ,C  ,G  ,A  ,T  ,A  ,
        /*e*/N  ,L  ,l  ,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/U  ,DUL,L  ,l  ,INF,INF,INF,INF,INF,INF,INF,
        /*C*/u  ,UL ,DUL,DUL,D  ,INF,INF,INF,INF,INF,INF,
        /*G*/INF,u  ,DUL,DUl,DUL,L  ,INF,INF,INF,INF,INF,
        /*A*/INF,INF,Du ,uL ,Ul ,DUL,L  ,INF,INF,INF,INF,
        /*C*/INF,INF,INF,Du ,uL ,Ul ,DUL,L  ,INF,INF,INF,
        /*T*/INF,INF,INF,INF,u  ,uL ,UL ,DUL,D  ,INF,INF,
        /*A*/INF,INF,INF,INF,INF,Du ,uL ,DUL,DUl,L  ,INF,
        /*G*/INF,INF,INF,INF,INF,INF,u  ,DuL,UL ,DUl,D  ,
        /*C*/INF,INF,INF,INF,INF,INF,INF,Du ,DuL,DUL,DUl
        }
    };
}();

static auto dna4_single_diagonal = []()
{
    //   012345|
    // 0 x     |
    // 1  x    |
    // 2   x   |
    // 3    x  |
    // 4     x |
    // 5      x|

    return alignment_fixture
    {
        "ATCGA"_dna4,
        "ACGAC"_dna4,
        align_config | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{0},
                                                          seqan3::align_cfg::upper_diagonal{0}} |
                       seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                    seqan3::mismatch_score{-5}}},
        -16,
        "ATCGA",
        "ACGAC",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 5u,
        /*.sequence2_end_position = */ 5u,
        std::vector<std::optional<int32_t>>
        {
        //     e,  A,  T,  C,  G,  A,
        /*e*/  0,INF,INF,INF,INF,INF,
        /*A*/INF,  4,INF,INF,INF,INF,
        /*C*/INF,INF, -1,INF,INF,INF,
        /*G*/INF,INF,INF, -6,INF,INF,
        /*A*/INF,INF,INF,INF,-11,INF,
        /*C*/INF,INF,INF,INF,INF,-16
        },
        std::vector<std::optional<seqan3::detail::trace_directions>>
        {
        //      e,  A,  T,  C,  G,  A,
        /*e*/   N,INF,INF,INF,INF,INF,
        /*A*/ INF,  D,INF,INF,INF,INF,
        /*C*/ INF,INF,  D,INF,INF,INF,
        /*G*/ INF,INF,INF,  D,INF,INF,
        /*A*/ INF,INF,INF,INF,  D,INF,
        /*C*/ INF,INF,INF,INF,INF,  D
        }
    };
}();

static auto dna4_large_band = []()
{
    //   012345|
    // 0       |
    // 1       |
    // 2       |
    // 3       |
    // 4       |
    // 5       |

    return alignment_fixture
    {
        "ATCGA"_dna4,
        "ACGAC"_dna4,
        align_config | seqan3::align_cfg::band_fixed_size{seqan3::align_cfg::lower_diagonal{-10},
                                                          seqan3::align_cfg::upper_diagonal{10}} |
                       seqan3::align_cfg::scoring_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                    seqan3::mismatch_score{-5}}},
        -6,
        "ATCGA-",
        "A-CGAC",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 5u,
        /*.sequence2_end_position = */ 5u,
        std::vector<std::optional<int32_t>>
        {
        //   e   ,A  ,T  ,C  ,G  ,A  ,
        /*e*/0  ,-11,-12,-13,-14,-15,
        /*A*/-11,4  ,-7 ,-8 ,-9 ,-10,
        /*C*/-12,-7 ,-1 ,-3 ,-13,-14,
        /*G*/-13,-8 ,-12,-6 ,1  ,-10,
        /*A*/-14,-9 ,-13,-15,-10,5  ,
        /*C*/-15,-10,-14,-9 ,-11,-6
        },
        std::vector<std::optional<seqan3::detail::trace_directions>>
        {
        //   e  ,A  ,T  ,C  ,G  ,A  ,
        /*e*/N  ,L  ,l  ,l  ,l  ,l  ,
        /*A*/U  ,DUL,L  ,l  ,l  ,DUl,
        /*C*/u  ,UL ,DUL,DUL,DUl,DUl,
        /*G*/u  ,uL ,DUL,DUl,DuL,L  ,
        /*A*/u  ,DuL,DuL,ul ,Ul ,DUL,
        /*C*/u  ,uL ,DuL,Dul,uL ,Ul
        }
    };
}();
// clang-format on

} // namespace seqan3::test::alignment::fixture::global::affine::banded
