// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <vector>

#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_mode.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "alignment_fixture.hpp"

using seqan3::operator""_dna4;

namespace seqan3::test::alignment::fixture::semi_global::affine::banded
{

inline constexpr auto align_config = seqan3::align_cfg::mode{seqan3::global_alignment} |
                                     seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1},
                                                                               seqan3::gap_open_score{-10}}} |
                                     seqan3::align_cfg::band{seqan3::static_band{seqan3::lower_bound{-4},
                                                                                 seqan3::upper_bound{8}}};

inline constexpr auto align_config_semi_seq1 = align_config | seqan3::align_cfg::aligned_ends{seqan3::free_ends_first};
inline constexpr auto align_config_semi_seq2 = align_config | seqan3::align_cfg::aligned_ends{seqan3::free_ends_second};

static auto dna4_01_semi_first = []()
{
    return alignment_fixture
    {
        "TTTTTACGTATGTCCCCC"_dna4,
        "ACGTAAAACGT"_dna4,
        align_config_semi_seq1
            | seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                           seqan3::mismatch_score{-5}}},
        10,
        "ACGT---ATGT",
        "ACGTAAAACGT",
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{5u}, seqan3::detail::row_index_type{0u}},
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{13u}, seqan3::detail::row_index_type{11u}},
        std::vector<std::optional<int32_t>>
        {
        //      e,  T,  T,  T,  T,  T,  A,  C,  G,  T,  A,  T,  G,  T,  C,  C,  C,  C,  C
        /*e*/  0 ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ -11,-5 ,-5 ,-5 ,-5 ,-5 ,4  ,-5 ,-5 ,-5 ,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ -12,-12,-10,-10,-10,-10,-7 ,8  ,-3 ,-4 ,-5 ,INF,INF,INF,INF,INF,INF,INF,INF,
        /*G*/ -13,-13,-13,-13,-13,-13,-8 ,-3 ,12 ,1  ,0  ,-1 ,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ -14,-9 ,-9 ,-9 ,-9 ,-9 ,-9 ,-4 ,1  ,16 ,5  ,4  ,3  ,INF,INF,INF,INF,INF,INF,
        /*A*/ INF,-15,-14,-14,-14,-14,-5 ,-5 ,0  ,5  ,20 ,9  ,8  ,7  ,INF,INF,INF,INF,INF,
        /*A*/ INF,INF,-16,-16,-16,-16,-10,-6 ,-1 ,4  ,9  ,15 ,4  ,3  ,2  ,INF,INF,INF,INF,
        /*A*/ INF,INF,INF,-17,-17,-17,-12,-7 ,-2 ,3  ,8  ,4  ,10 ,-1 ,-2 ,-3 ,INF,INF,INF,
        /*A*/ INF,INF,INF,INF,-18,-18,-13,-8 ,-3 ,2  ,7  ,3  ,-1 ,5  ,-6 ,-7 ,-8 ,INF,INF,
        /*C*/ INF,INF,INF,INF,INF,-19,-14,-9 ,-4 ,1  ,6  ,2  ,-2 ,-6 ,9  ,-2 ,-3 ,-4 ,INF,
        /*G*/ INF,INF,INF,INF,INF,INF,-15,-10,-5 ,0  ,5  ,1  ,6  ,-5 ,-2 ,4  ,-7 ,-8 ,-9 ,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,-11,-6 ,-1 ,4  ,9  ,-2 ,10 ,-1 ,-2 ,-1 ,-4 ,-5
        },
        std::vector<std::optional<seqan3::detail::trace_directions>>
        {
        //      e,  T,  T,  T,  T,  T,  A,  C,  G,  T,  A,  T,  G,  T,  C,  C,  C,  C,  C
        /*e*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ U  ,DUL,DUL,DUL,DUL,DUL,DUL,DUL,DUl,D  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ u  ,uL ,DuL,DuL,DuL,DuL,UL ,DuL,L  ,l  ,l  ,INF,INF,INF,INF,INF,INF,INF,INF,
        /*G*/ u  ,uL ,uL ,uL ,uL ,uL ,uL ,UL ,DuL,L  ,l  ,l  ,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ u  ,DuL,DuL,DuL,DuL,DuL,uL ,uL ,UL ,DUL,L  ,DUl,l  ,INF,INF,INF,INF,INF,INF,
        /*A*/ INF,u  ,DuL,DuL,DuL,DuL,DuL,uL ,uL ,UL ,DUL,L  ,l  ,l  ,INF,INF,INF,INF,INF,
        /*A*/ INF,INF,u  ,uL ,uL ,uL ,DuL,uL ,uL ,uL ,DUL,DUL,DUL,DUl,D  ,INF,INF,INF,INF,
        /*A*/ INF,INF,INF,u  ,uL ,uL ,DuL,uL ,uL ,uL ,DuL,DUL,Dul,DuL,DUl,D  ,INF,INF,INF,
        /*A*/ INF,INF,INF,INF,u  ,uL ,DuL,uL ,uL ,uL ,DuL,DuL,DUl,Dul,DuL,DUl,D  ,INF,INF,
        /*C*/ INF,INF,INF,INF,INF,u  ,uL ,DuL,uL ,uL ,uL ,DuL,Dul,DUl,Dul,DuL,DUl,D  ,INF,
        /*G*/ INF,INF,INF,INF,INF,INF,u  ,uL ,DuL,uL ,uL ,DuL,Dul,L  ,Ul ,DUl,DUL,DUl,D  ,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,u  ,uL ,DuL,uL ,DuL,L  ,Dul,L  ,l  ,Dul,l  ,l
        }
    };
}();

// TODO: Fix the computation of the trailing gap, which depends on whether the score in the last cell
// was a gap or not (the gap must be extended in this case.) -> or throw invalid_alignment_configuration if
// the score cannot be computed.
// static auto dna4_02_semi_first = []()
// {
//     return alignment_fixture
//     {
//         "ACGTAAAACGT"_dna4,
//         "TTTTTACGTATGTCCCCC"_dna4,
//         align_config_semi_seq1
//             | seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
//                                                                            seqan3::mismatch_score{-5}}},
//         -13,
//         "-----ACGTA--------AAACGT",
//         "TTTTTACGTATGTCCCCC------",
//         seqan3::alignment_coordinate{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}},
//         seqan3::alignment_coordinate{seqan3::detail::column_index_type{5u}, seqan3::detail::row_index_type{18u}}
//     };
// }();

static auto dna4_03_semi_second = []()
{
    return alignment_fixture
    {
        "TTTTTACGTATGTCCCCC"_dna4,
        "ACGTAAAACGT"_dna4,
        align_config_semi_seq2
            | seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                           seqan3::mismatch_score{-5}}},
        -19,
        "TTTTTACGTATGTCCCCC",
        "GTAAAACGT---------",
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{2u}},
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{18u}, seqan3::detail::row_index_type{11u}},
        std::vector<std::optional<int32_t>>
        {
        //      e,  T,  T,  T,  T,  T,  A,  C,  G,  T,  A,  T,  G,  T,  C,  C,  C,  C,  C
        /*e*/ 0  ,-11,-12,-13,-14,-15,-16,-17,-18,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ 0  ,-5 ,-12,-13,-14,-15,-11,-17,-18,-19,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ 0  ,-5 ,-10,-13,-14,-15,-16,-7 ,-18,-19,-20,INF,INF,INF,INF,INF,INF,INF,INF,
        /*G*/ 0  ,-5 ,-10,-13,-14,-15,-16,-17,-3 ,-14,-15,-16,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ 0  ,4  ,-1 ,-6 ,-9 ,-10,-11,-12,-13,1  ,-10,-11,-12,INF,INF,INF,INF,INF,INF,
        /*A*/ INF,-5 ,-1 ,-6 ,-11,-14,-6 ,-16,-15,-10,5  ,-6 ,-7 ,-8 ,INF,INF,INF,INF,INF,
        /*A*/ INF,INF,-10,-6 ,-11,-16,-10,-11,-16,-11,-6 ,0  ,-11,-12,-13,INF,INF,INF,INF,
        /*A*/ INF,INF,INF,-15,-11,-16,-12,-15,-16,-12,-7 ,-11,-5 ,-16,-17,-18,INF,INF,INF,
        /*A*/ INF,INF,INF,INF,-20,-16,-12,-17,-18,-13,-8 ,-12,-16,-10,-21,-22,-23,INF,INF,
        /*C*/ INF,INF,INF,INF,INF,-25,-20,-8 ,-19,-14,-9 ,-13,-17,-21,-6 ,-17,-18,-19,INF,
        /*G*/ INF,INF,INF,INF,INF,INF,-21,-19,-4 ,-15,-10,-14,-9 ,-19,-17,-11,-22,-23,-24,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,-20,-15,0  ,-11,-6 ,-13,-5 ,-15,-16,-16,-18,-19
        },
        std::vector<std::optional<seqan3::detail::trace_directions>>
        {
        //      e,  T,  T,  T,  T,  T,  A,  C,  G,  T,  A,  T,  G,  T,  C,  C,  C,  C,  C
        /*e*/ N  ,L  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ N  ,DUL,l  ,l  ,l  ,l  ,DUl,l  ,l  ,l  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ N  ,DUL,DUl,l  ,l  ,l  ,l  ,DUl,l  ,l  ,l  ,INF,INF,INF,INF,INF,INF,INF,INF,
        /*G*/ N  ,DUL,DUl,l  ,l  ,l  ,l  ,l  ,DUl,L  ,l  ,l  ,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ N  ,DUL,DUL,DUl,DUl,DUl,l  ,l  ,l  ,DUl,L  ,DUl,l  ,INF,INF,INF,INF,INF,INF,
        /*A*/ INF,DU ,DUL,DUL,DUl,DUl,DUl,Dul,ul ,Ul ,DUl,L  ,l  ,l  ,INF,INF,INF,INF,INF,
        /*A*/ INF,INF,DU ,DUL,DuL,Dul,DUl,Dul,ul ,ul ,DUL,DUL,DUL,DUl,D  ,INF,INF,INF,INF,
        /*A*/ INF,INF,INF,DU ,DuL,DuL,Dul,DuL,Dul,ul ,DuL,DUL,Dul,DuL,DUl,D  ,INF,INF,INF,
        /*A*/ INF,INF,INF,INF,DU ,DuL,DuL,DuL,ul ,ul ,DuL,DuL,DUl,Dul,DuL,DUl,D  ,INF,INF,
        /*C*/ INF,INF,INF,INF,INF,Du ,uL ,DuL,uL ,ul ,ul ,DuL,Dul,DUl,Dul,DuL,DUl,D  ,INF,
        /*G*/ INF,INF,INF,INF,INF,INF,u  ,UL ,DuL,uL ,ul ,Dul,Dul,l  ,Ul ,DUl,DUl,DUl,D  ,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,u  ,UL ,DuL,uL ,Dul,l  ,Dul,l  ,l  ,Dul,l  ,l
        }
    };
}();

static auto dna4_04_semi_second = []()
{
    return alignment_fixture
    {
        "ACGTAAAACGT"_dna4,
        "TTTTTACGTATGTCCCCC"_dna4,
        align_config_semi_seq2
            | seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                           seqan3::mismatch_score{-5}}},
        -5,
        "ACGTAAAACGT",
        "------TACGT",
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{4u}},
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{11u}, seqan3::detail::row_index_type{9u}},
        std::vector<std::optional<int32_t>>
        {
        //      e,  A,  C,  G,  T,  A,  A,  A,  A,  C,  G,  T
        /*e*/ 0  ,-11,-12,-13,-14,-15,-16,-17,-18,INF,INF,INF,
        /*T*/ 0  ,-5 ,-12,-13,-9 ,-15,-16,-17,-18,-19,INF,INF,
        /*T*/ 0  ,-5 ,-10,-13,-9 ,-14,-16,-17,-18,-19,-20,INF,
        /*T*/ 0  ,-5 ,-10,-13,-9 ,-14,-16,-17,-18,-19,-20,-16,
        /*T*/ 0  ,-5 ,-10,-13,-9 ,-14,-16,-17,-18,-19,-20,-16,
        /*T*/ INF,-5 ,-10,-15,-9 ,-14,-19,-21,-22,-23,-24,-16,
        /*A*/ INF,INF,-10,-15,-20,-5 ,-10,-15,-17,-19,-20,-21,
        /*C*/ INF,INF,INF,-15,-20,-16,-10,-15,-20,-13,-24,-25,
        /*G*/ INF,INF,INF,INF,-20,-17,-21,-15,-20,-24,-9 ,-20,
        /*T*/ INF,INF,INF,INF,INF,-18,-22,-26,-20,-25,-20,-5 ,
        /*A*/ INF,INF,INF,INF,INF,INF,-14,-18,-22,-25,-21,-16,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,-19,-23,-27,-22,-17,
        /*G*/ INF,INF,INF,INF,INF,INF,INF,INF,-24,-28,-23,-18,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,-29,-24,-19,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,-25,-20,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,-21,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF
        },
        std::vector<std::optional<seqan3::detail::trace_directions>>
        {
        //      e,  A,  C,  G,  T,  A,  A,  A,  A,  C,  G,  T
        /*e*/ N  ,L  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,INF,INF,INF,
        /*T*/ N  ,DUL,l  ,l  ,DUl,l  ,l  ,l  ,l  ,l  ,INF,INF,
        /*T*/ N  ,DUL,DUl,l  ,DUl,DUl,l  ,l  ,l  ,l  ,l  ,INF,
        /*T*/ N  ,DUL,DUl,l  ,DUl,DUl,l  ,l  ,l  ,l  ,l  ,D  ,
        /*T*/ N  ,DUL,DUl,l  ,DUl,DUl,l  ,l  ,l  ,l  ,l  ,DUl,
        /*T*/ INF,DU ,DUL,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,
        /*A*/ INF,INF,DU ,DuL,DUl,DUl,DuL,Dul,Dul,l  ,l  ,l  ,
        /*C*/ INF,INF,INF,Du ,DuL,Ul ,DUL,DUL,DUl,DUl,DUl,Dul,
        /*G*/ INF,INF,INF,INF,Du ,uL ,DUL,DUl,DuL,Ul ,Dul,L  ,
        /*T*/ INF,INF,INF,INF,INF,u  ,DuL,DUl,Dul,DuL,Ul ,DuL,
        /*A*/ INF,INF,INF,INF,INF,INF,Du ,DuL,Dul,Dul,ul ,Ul ,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,Du ,DuL,Dul,ul ,Dul,
        /*G*/ INF,INF,INF,INF,INF,INF,INF,INF,Du ,DuL,Dul,uL ,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,Du ,uL ,DuL,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,u  ,uL ,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,u  ,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF
        }
    };
}();

} // namespace seqan3::test::alignment::fixture::global::affine::banded
