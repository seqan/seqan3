// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <vector>

#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_mode.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "alignment_fixture.hpp"

namespace seqan3::test::alignment::fixture::global::affine::unbanded
{

using namespace seqan3;
using namespace seqan3::detail;

inline constexpr auto align_config = align_cfg::mode{global_alignment} |
                                     align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}};

static auto dna4_01 = []()
{
    return alignment_fixture
    {
        // score: 8 (7 insertions, 1 substitutions)
        // alignment:
        // AACCGGTTAACCGGTT
        // | | | | | | | |
        // A-C-G-T-A-C-G-TA
        "AACCGGTTAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        align_config | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        -18,
        "AACCGGTTAACCG---GTT",
        "A----------CGTACGTA",
        alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{16u}, row_index_type{9u}},
        std::vector
        {
         //        A  ,A  ,C  ,C  ,G  ,G  ,T  ,T  ,A  ,A  ,C  ,C  ,G  ,G  ,T  ,T
                0 ,-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,-22,-23,-24,-25,-26,
        /*A*/  -11,  4, -7, -8, -9,-10,-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,
        /*C*/  -12, -7, -1, -3, -4,-14,-15,-16,-17,-18,-19,-11,-12,-22,-23,-24,-25,
        /*G*/  -13, -8,-12, -6, -8,  0,-10,-12,-13,-14,-15,-16,-16,-8 ,-18,-20,-21,
        /*T*/  -14, -9,-13,-15,-11,-11,-5 , -6, -8,-18,-19,-20,-21,-19,-13,-14,-16,
        /*A*/  -15,-10,-5 ,-16,-17,-12,-16,-10,-11, -4,-14,-16,-17,-18,-19,-18,-19,
        /*C*/  -16,-11,-15, -1,-12,-13,-14,-15,-15,-15, -9,-10,-12,-21,-22,-23,-23,
        /*G*/  -17,-12,-16,-12, -6, -8, -9,-19,-20,-16,-20,-14,-15, -8,-17,-20,-21,
        /*T*/  -18,-13,-17,-13,-17,-11,-13, -5,-15,-17,-18,-19,-19,-19,-13,-13,-16,
        /*A*/  -19,-14,-9 ,-14,-18,-16,-16,-16,-10,-11,-13,-23,-24,-20,-24,-18,-18
        },
        std::vector
        {
        //           A  ,A  ,C  ,C  ,G  ,G  ,T  ,T  ,A  ,A  ,C  ,C  ,G  ,G  ,T  ,T
        /*e*/    N  ,L  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*A*/    U  ,DUL,DUL,l  ,l  ,l  ,l  ,l  ,l  ,DUl,DUl,l  ,l  ,l  ,l  ,l  ,l  ,
        /*C*/    u  ,UL ,DUL,DUL,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,
        /*G*/    u  ,uL ,DUL,DUl,DUL,Dul,DuL,l  ,l  ,l  ,l  ,l  ,DUl,Dul,Dul,l  ,l  ,
        /*T*/    u  ,uL ,DuL,ul ,Dul,UL ,DUL,DUL,DUl,DUl,DUl,Dul,Dul,Ul ,DUl,DUl,DUl,
        /*A*/    u  ,DuL,DuL,uL ,ul ,ul ,DUl,DUl,DUl,Dul,DuL,l  ,l  ,l  ,l  ,DUl,DUl,
        /*C*/    u  ,uL ,DuL,Dul,DuL,ul ,l  ,l  ,Dul,Ul ,DUl,Dul,Dul,ul ,l  ,l  ,Dul,
        /*G*/    u  ,uL ,DuL,Ul ,DuL,DuL,Dul,Dul,Dul,ul ,DUl,DUl,DUl,Dul,DuL,l  ,l  ,
        /*T*/    u  ,uL ,DuL,ul ,DUL,Dul,DuL,Dul,DuL,ul ,l  ,l  ,Dul,Ul ,Dul,Dul,Dul,
        /*A*/    u  ,DuL,DuL,uL ,Dul,ul ,Dul,Ul ,Dul,DuL,Dul,Dul,Dul,ul ,DUl,DUl,DUl
        }
    };
}();

static auto dna4_02 = []()
{
    return alignment_fixture
    {
        "ACGTACGTA"_dna4,
        "AACCGGTTAACCGGTT"_dna4,
        align_config | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        -18,
        "ACGTAC----------GTA",
        "A---ACCGGTTAACCGGTT",
        alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{9u}, row_index_type{16u}},
        std::vector
        {
        //		      A,  C,  G,  T,  A,  C,  G,  T,  A
                  0,-11,-12,-13,-14,-15,-16,-17,-18,-19,
        /*A*/	-11,  4, -7, -8, -9,-10,-11,-12,-13,-14,
        /*A*/	-12, -7, -1,-12,-13, -5,-15,-16,-17, -9,
        /*C*/	-13, -8, -3, -6,-15,-16, -1,-12,-13,-14,
        /*C*/	-14, -9, -4, -8,-11,-17,-12, -6,-17,-18,
        /*G*/	-15,-10,-14,  0,-11,-12,-13, -8,-11,-16,
        /*G*/	-16,-11,-15,-10, -5,-16,-14, -9,-13,-16,
        /*T*/	-17,-12,-16,-12, -6,-10,-15,-19, -5,-16,
        /*T*/	-18,-13,-17,-13, -8,-11,-15,-20,-15,-10,
        /*A*/	-19,-14,-18,-14,-18, -4,-15,-16,-17,-11,
        /*A*/	-20,-15,-19,-15,-19,-14, -9,-20,-18,-13,
        /*C*/	-21,-16,-11,-16,-20,-16,-10,-14,-19,-23,
        /*C*/	-22,-17,-12,-16,-21,-17,-12,-15,-19,-24,
        /*G*/	-23,-18,-22, -8,-19,-18,-21, -8,-19,-20,
        /*G*/	-24,-19,-23,-18,-13,-19,-22,-17,-13,-24,
        /*T*/	-25,-20,-24,-20,-14,-18,-23,-20,-13,-18,
        /*T*/	-26,-21,-25,-21,-16,-19,-23,-21,-16,-18
        },
        std::vector
        {
        //		      A,  C,  G,  T,  A,  C,  G,  T,  A
        /*e*/   N  ,L  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*A*/   U  ,DUL,L  ,l  ,l  ,DUl,l  ,l  ,l  ,DUl,
        /*A*/   u  ,DUL,DUL,DUL,DUl,DUl,DUl,DUl,DUl,DUl,
        /*C*/   u  ,uL ,DUL,DuL,l  ,Ul ,Dul,L  ,l  ,l  ,
        /*C*/   u  ,uL ,DuL,DUL,Dul,ul ,DUl,DUl,DUL,Dul,
        /*G*/   u  ,uL ,DuL,Dul,L  ,l  ,ul ,DUl,Dul,l  ,
        /*G*/   u  ,uL ,DuL,DUl,DUL,DuL,ul ,Dul,DUl,Dul,
        /*T*/   u  ,uL ,DuL,ul ,DUL,DuL,ul ,Dul,Dul,L  ,
        /*T*/   u  ,uL ,DuL,ul ,DuL,DuL,Dul,Dul,DUl,Dul,
        /*A*/   u  ,DuL,DuL,ul ,DuL,Dul,L  ,l  ,ul ,DUl,
        /*A*/   u  ,DuL,DuL,ul ,DuL,DUl,DuL,DuL,ul ,Dul,
        /*C*/   u  ,uL ,DuL,uL ,Dul,ul ,Dul,DuL,ul ,Dul,
        /*C*/   u  ,uL ,DuL,DuL,Dul,ul ,Dul,DuL,Dul,Dul,
        /*G*/   u  ,uL ,DuL,Dul,L  ,ul ,ul ,Dul,L  ,l  ,
        /*G*/   u  ,uL ,DuL,Dul,DuL,uL ,ul ,DUl,Dul,DuL,
        /*T*/   u  ,uL ,DuL,ul ,DuL,DuL,ul ,ul ,Dul,DuL,
        /*T*/   u  ,uL ,DuL,ul ,DuL,DuL,Dul,ul ,Dul,DuL,



















        }
    };
}();

} // namespace seqan3::test::alignment::fixture::global::affine::unbanded
