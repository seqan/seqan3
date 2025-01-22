// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

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

namespace seqan3::test::alignment::fixture::global::affine::unbanded
{
// clang-format off

inline constexpr auto align_config = seqan3::align_cfg::method_global{}
                                   | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-10},
                                                                        seqan3::align_cfg::extension_score{-1}};
inline constexpr auto align_config_dna_score = align_config
                                             | seqan3::align_cfg::scoring_scheme{
                                                    seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                      seqan3::mismatch_score{-5}}};

// Naming is used to disambiguate configs for global alignment with different scoring parameters.
// `part_XX` indicates that the config is the same but the sequences differ in their lengths.
// This is helpful to collect different test fixtures to test vectorised alignments where we can
// only apply alignments with the same config.
static auto dna4_match_4_mismatch_5_gap_1_open_10_part_01 = []()
{
    using seqan3::operator""_dna4;
    return alignment_fixture
    {
        // score: 8 (7 insertions, 1 substitutions)
        // alignment:
        // AACCGGTTAACCGGTT
        // | | | | | | | |
        // A-C-G-T-A-C-G-TA
        "AACCGGTTAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        align_config_dna_score,
        -18,
        "AACCGGTTAACCG---GTT",
        "A----------CGTACGTA",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 16u,
        /*.sequence2_end_position = */ 9u,
        std::vector
        {
        //      e,  A,  A,  C,  C,  G,  G,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/  0 ,-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,-22,-23,-24,-25,-26,
        /*A*/ -11,  4, -7, -8, -9,-10,-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,
        /*C*/ -12, -7, -1, -3, -4,-14,-15,-16,-17,-18,-19,-11,-12,-22,-23,-24,-25,
        /*G*/ -13, -8,-12, -6, -8,  0,-10,-12,-13,-14,-15,-16,-16,-8 ,-18,-20,-21,
        /*T*/ -14, -9,-13,-15,-11,-11,-5 , -6, -8,-18,-19,-20,-21,-19,-13,-14,-16,
        /*A*/ -15,-10,-5 ,-16,-17,-12,-16,-10,-11, -4,-14,-16,-17,-18,-19,-18,-19,
        /*C*/ -16,-11,-15, -1,-12,-13,-14,-15,-15,-15, -9,-10,-12,-21,-22,-23,-23,
        /*G*/ -17,-12,-16,-12, -6, -8, -9,-19,-20,-16,-20,-14,-15, -8,-17,-20,-21,
        /*T*/ -18,-13,-17,-13,-17,-11,-13, -5,-15,-17,-18,-19,-19,-19,-13,-13,-16,
        /*A*/ -19,-14,-9 ,-14,-18,-16,-16,-16,-10,-11,-13,-23,-24,-20,-24,-18,-18
        },
        std::vector
        {
        //      e,  A,  A,  C,  C,  G,  G,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/ N  ,L  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*A*/ U  ,DUL,DUL,l  ,l  ,l  ,l  ,l  ,l  ,DUl,DUl,l  ,l  ,l  ,l  ,l  ,l  ,
        /*C*/ u  ,UL ,DUL,DUL,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,
        /*G*/ u  ,uL ,DUL,DUl,DUL,Dul,DuL,l  ,l  ,l  ,l  ,l  ,DUl,Dul,Dul,l  ,l  ,
        /*T*/ u  ,uL ,DuL,ul ,Dul,UL ,DUL,DUL,DUl,DUl,DUl,Dul,Dul,Ul ,DUl,DUl,DUl,
        /*A*/ u  ,DuL,DuL,uL ,ul ,ul ,DUl,DUl,DUl,Dul,DuL,l  ,l  ,l  ,l  ,DUl,DUl,
        /*C*/ u  ,uL ,DuL,Dul,DuL,ul ,l  ,l  ,Dul,Ul ,DUl,Dul,Dul,ul ,l  ,l  ,Dul,
        /*G*/ u  ,uL ,DuL,Ul ,DuL,DuL,Dul,Dul,Dul,ul ,DUl,DUl,DUl,Dul,DuL,l  ,l  ,
        /*T*/ u  ,uL ,DuL,ul ,DUL,Dul,DuL,Dul,DuL,ul ,l  ,l  ,Dul,Ul ,Dul,Dul,Dul,
        /*A*/ u  ,DuL,DuL,uL ,Dul,ul ,Dul,Ul ,Dul,DuL,Dul,Dul,Dul,ul ,DUl,DUl,DUl
        }
    };
}();

static auto dna4_match_4_mismatch_5_gap_1_open_10_part_02 = []()
{
    using seqan3::operator""_dna4;
    return alignment_fixture
    {
        "ACGTACGTA"_dna4,
        "AACCGGTTAACCGGTT"_dna4,
        align_config_dna_score,
        -18,
        "ACGTAC----------GTA",
        "A---ACCGGTTAACCGGTT",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 9u,
        /*.sequence2_end_position = */ 16u,
        std::vector
        {
        //      e,  A,  C,  G,  T,  A,  C,  G,  T,  A
        /*e*/ 0  ,-11,-12,-13,-14,-15,-16,-17,-18,-19,
        /*A*/ -11,  4, -7, -8, -9,-10,-11,-12,-13,-14,
        /*A*/ -12, -7, -1,-12,-13, -5,-15,-16,-17, -9,
        /*C*/ -13, -8, -3, -6,-15,-16, -1,-12,-13,-14,
        /*C*/ -14, -9, -4, -8,-11,-17,-12, -6,-17,-18,
        /*G*/ -15,-10,-14,  0,-11,-12,-13, -8,-11,-16,
        /*G*/ -16,-11,-15,-10, -5,-16,-14, -9,-13,-16,
        /*T*/ -17,-12,-16,-12, -6,-10,-15,-19, -5,-16,
        /*T*/ -18,-13,-17,-13, -8,-11,-15,-20,-15,-10,
        /*A*/ -19,-14,-18,-14,-18, -4,-15,-16,-17,-11,
        /*A*/ -20,-15,-19,-15,-19,-14, -9,-20,-18,-13,
        /*C*/ -21,-16,-11,-16,-20,-16,-10,-14,-19,-23,
        /*C*/ -22,-17,-12,-16,-21,-17,-12,-15,-19,-24,
        /*G*/ -23,-18,-22, -8,-19,-18,-21, -8,-19,-20,
        /*G*/ -24,-19,-23,-18,-13,-19,-22,-17,-13,-24,
        /*T*/ -25,-20,-24,-20,-14,-18,-23,-20,-13,-18,
        /*T*/ -26,-21,-25,-21,-16,-19,-23,-21,-16,-18
        },
        std::vector
        {
        //	    e,  A,  C,  G,  T,  A,  C,  G,  T,  A
        /*e*/ N  ,L  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*A*/ U  ,DUL,L  ,l  ,l  ,DUl,l  ,l  ,l  ,DUl,
        /*A*/ u  ,DUL,DUL,DUL,DUl,DUl,DUl,DUl,DUl,DUl,
        /*C*/ u  ,uL ,DUL,DuL,l  ,Ul ,Dul,L  ,l  ,l  ,
        /*C*/ u  ,uL ,DuL,DUL,Dul,ul ,DUl,DUl,DUL,Dul,
        /*G*/ u  ,uL ,DuL,Dul,L  ,l  ,ul ,DUl,Dul,l  ,
        /*G*/ u  ,uL ,DuL,DUl,DUL,DuL,ul ,Dul,DUl,Dul,
        /*T*/ u  ,uL ,DuL,ul ,DUL,DuL,ul ,Dul,Dul,L  ,
        /*T*/ u  ,uL ,DuL,ul ,DuL,DuL,Dul,Dul,DUl,Dul,
        /*A*/ u  ,DuL,DuL,ul ,DuL,Dul,L  ,l  ,ul ,DUl,
        /*A*/ u  ,DuL,DuL,ul ,DuL,DUl,DuL,DuL,ul ,Dul,
        /*C*/ u  ,uL ,DuL,uL ,Dul,ul ,Dul,DuL,ul ,Dul,
        /*C*/ u  ,uL ,DuL,DuL,Dul,ul ,Dul,DuL,Dul,Dul,
        /*G*/ u  ,uL ,DuL,Dul,L  ,ul ,ul ,Dul,L  ,l  ,
        /*G*/ u  ,uL ,DuL,Dul,DuL,uL ,ul ,DUl,Dul,DuL,
        /*T*/ u  ,uL ,DuL,ul ,DuL,DuL,ul ,ul ,Dul,DuL,
        /*T*/ u  ,uL ,DuL,ul ,DuL,DuL,Dul,ul ,Dul,DuL,
        }
    };
}();

static auto dna4_match_4_mismatch_5_gap_1_open_10_part_03 = []()
{
    using seqan3::operator""_dna4;
    return alignment_fixture
    {
        "AGATCGACTAGCGAGCTACGAGCTAGC"_dna4,
        "AGACGATCGACGAGCGACTACGTACGA"_dna4,
        align_config_dna_score,
        26,
        "A---GATCGACTAGCGAGCTACGAGCTA-GC",
        "AGACGATCGACGAGCGA-CTACG---TACGA",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 27u,
        /*.sequence2_end_position = */ 27u,
        std::vector
        {
        //      e,  A,  G,  A,  T,  C,  G,  A,  C,  T,  A,  G,  C,  G,  A,  G,  C,  T,  A,  C,  G,  A,  G,  C,  T,  A,  G,  C,
        /*e*/ 0  ,-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,-22,-23,-24,-25,-26,-27,-28,-29,-30,-31,-32,-33,-34,-35,-36,-37,
        /*A*/ -11,4  ,-7 ,-8 ,-9 ,-10,-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,-22,-23,-24,-25,-26,-27,-28,-29,-30,-31,-32,
        /*G*/ -12,-7 ,8  ,-3 ,-4 ,-5 ,-6 ,-7 ,-8 ,-9 ,-10,-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,-22,-23,-24,-25,-26,-27,
        /*A*/ -13,-8 ,-3 ,12 ,1  ,0  ,-1 ,-2 ,-3 ,-4 ,-5 ,-6 ,-7 ,-8 ,-9 ,-10,-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,-22,
        /*C*/ -14,-9 ,-4 ,1  ,7  ,5  ,-5 ,-6 ,2  ,-8 ,-9 ,-10,-2 ,-12,-13,-14,-6 ,-16,-17,-9 ,-19,-20,-21,-13,-23,-24,-25,-17,
        /*G*/ -15,-10,-5 ,0  ,-4 ,2  ,9  ,-2 ,-3 ,-3 ,-5 ,-5 ,-7 ,2  ,-9 ,-9 ,-11,-11,-13,-14,-5 ,-16,-16,-18,-18,-20,-20,-22,
        /*A*/ -16,-11,-6 ,-1 ,-5 ,-7 ,-2 ,13 ,2  ,1  ,1  ,-1 ,-2 ,-3 ,6  ,-5 ,-6 ,-7 ,-7 ,-9 ,-10,-1 ,-12,-13,-14,-14,-16,-17,
        /*T*/ -17,-12,-7 ,-2 ,3  ,-8 ,-3 ,2  ,8  ,6  ,-4 ,-4 ,-6 ,-7 ,-5 ,1  ,-10,-2 ,-12,-12,-14,-12,-6 ,-17,-9 ,-19,-19,-21,
        /*C*/ -18,-13,-8 ,-3 ,-7 ,7  ,-4 ,1  ,6  ,3  ,1  ,-7 ,0  ,-9 ,-6 ,-10,5  ,-6 ,-7 ,-8 ,-9 ,-10,-11,-2 ,-13,-14,-15,-15,
        /*G*/ -19,-14,-9 ,-4 ,-8 ,-4 ,11 ,0  ,-1 ,1  ,-2 ,5  ,-5 ,4  ,-7 ,-2 ,-6 ,0  ,-11,-12,-4 ,-14,-6 ,-13,-7 ,-18,-10,-20,
        /*A*/ -20,-15,-10,-5 ,-9 ,-5 ,0  ,15 ,4  ,3  ,5  ,1  ,0  ,-1 ,8  ,-3 ,-4 ,-5 ,4  ,-7 ,-8 ,0  ,-10,-11,-12,-3 ,-14,-15,
        /*C*/ -21,-16,-11,-6 ,-10,-5 ,-1 ,4  ,19 ,8  ,7  ,6  ,5  ,4  ,3  ,3  ,1  ,0  ,-1 ,8  ,-3 ,-4 ,-5 ,-6 ,-7 ,-8 ,-8 ,-10,
        /*G*/ -22,-17,-12,-7 ,-11,-7 ,-1 ,3  ,8  ,14 ,3  ,11 ,1  ,9  ,-1 ,7  ,-2 ,-4 ,-5 ,-3 ,12 ,1  ,0  ,-1 ,-2 ,-3 ,-4 ,-5 ,
        /*A*/ -23,-18,-13,-8 ,-12,-8 ,-3 ,3  ,7  ,3  ,18 ,7  ,6  ,5  ,13 ,3  ,2  ,1  ,0  ,-1 ,1  ,16 ,5  ,4  ,3  ,2  ,1  ,0  ,
        /*G*/ -24,-19,-14,-9 ,-13,-9 ,-4 ,1  ,6  ,2  ,7  ,22 ,11 ,10 ,9  ,17 ,7  ,6  ,5  ,4  ,3  ,5  ,20 ,9  ,8  ,7  ,6  ,5  ,
        /*C*/ -25,-20,-15,-10,-14,-9 ,-5 ,0  ,5  ,1  ,6  ,11 ,26 ,15 ,14 ,13 ,21 ,11 ,10 ,9  ,8  ,7  ,9  ,24 ,13 ,12 ,11 ,10 ,
        /*G*/ -26,-21,-16,-11,-15,-11,-5 ,-1 ,4  ,0  ,5  ,10 ,15 ,30 ,19 ,18 ,17 ,16 ,15 ,14 ,13 ,12 ,11 ,13 ,19 ,8  ,16 ,6  ,
        /*A*/ -27,-22,-17,-12,-16,-12,-7 ,-1 ,3  ,-1 ,4  ,9  ,14 ,19 ,34 ,23 ,22 ,21 ,20 ,19 ,18 ,17 ,16 ,15 ,14 ,23 ,12 ,11 ,
        /*C*/ -28,-23,-18,-13,-17,-12,-8 ,-3 ,3  ,-2 ,3  ,8  ,13 ,18 ,23 ,29 ,27 ,17 ,16 ,24 ,14 ,13 ,12 ,20 ,10 ,12 ,18 ,16 ,
        /*T*/ -29,-24,-19,-14,-9 ,-14,-9 ,-4 ,1  ,7  ,2  ,7  ,12 ,17 ,22 ,18 ,24 ,31 ,20 ,19 ,19 ,17 ,16 ,15 ,24 ,13 ,12 ,13 ,
        /*A*/ -30,-25,-20,-15,-19,-14,-10,-5 ,0  ,-4 ,11 ,6  ,11 ,16 ,21 ,17 ,15 ,20 ,35 ,24 ,23 ,23 ,21 ,20 ,19 ,28 ,17 ,16 ,
        /*C*/ -31,-26,-21,-16,-20,-15,-11,-6 ,-1 ,-5 ,0  ,6  ,10 ,15 ,20 ,16 ,21 ,19 ,24 ,39 ,28 ,27 ,26 ,25 ,24 ,23 ,23 ,21 ,
        /*G*/ -32,-27,-22,-17,-21,-17,-11,-7 ,-2 ,-6 ,-1 ,4  ,9  ,14 ,19 ,24 ,13 ,18 ,23 ,28 ,43 ,32 ,31 ,30 ,29 ,28 ,27 ,26 ,
        /*T*/ -33,-28,-23,-18,-13,-18,-13,-8 ,-3 ,2  ,-2 ,3  ,8  ,13 ,18 ,14 ,19 ,17 ,22 ,27 ,32 ,38 ,27 ,26 ,34 ,24 ,23 ,22 ,
        /*A*/ -34,-29,-24,-19,-23,-18,-14,-9 ,-4 ,-8 ,6  ,2  ,7  ,12 ,17 ,13 ,11 ,16 ,21 ,26 ,31 ,36 ,33 ,24 ,23 ,38 ,27 ,26 ,
        /*C*/ -35,-30,-25,-20,-24,-19,-15,-10,-5 ,-9 ,-4 ,1  ,6  ,11 ,16 ,12 ,17 ,15 ,20 ,25 ,30 ,26 ,31 ,37 ,26 ,27 ,33 ,31 ,
        /*G*/ -36,-31,-26,-21,-25,-21,-15,-11,-6 ,-10,-5 ,0  ,5  ,10 ,15 ,20 ,9  ,14 ,19 ,24 ,29 ,25 ,30 ,26 ,32 ,26 ,31 ,28 ,
        /*A*/ -37,-32,-27,-22,-26,-22,-17,-11,-7 ,-11,-6 ,-1 ,4  ,9  ,14 ,10 ,15 ,13 ,18 ,23 ,28 ,33 ,22 ,25 ,21 ,36 ,25 ,26
        },
        std::vector
        {
        //      e,  A,  G,  A,  T,  C,  G,  A,  C,  T,  A,  G,  C,  G,  A,  G,  C,  T,  A,  C,  G,  A,  G,  C,  T,  A,  G,  C,
        /*e*/ N  ,L  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*A*/ U  ,DUL,L  ,DUl,l  ,l  ,l  ,DUl,l  ,l  ,DUl,l  ,l  ,l  ,DUl,l  ,l  ,l  ,DUl,l  ,l  ,DUl,l  ,l  ,l  ,DUl,l  ,l  ,
        /*G*/ u  ,UL ,DUL,L  ,l  ,l  ,DUl,l  ,l  ,l  ,l  ,DUl,l  ,DUl,l  ,DUl,l  ,l  ,l  ,l  ,DUl,l  ,DUl,l  ,l  ,l  ,DUl,l  ,
        /*A*/ u  ,DuL,UL ,DUL,L  ,l  ,l  ,DUl,l  ,l  ,DUl,l  ,l  ,l  ,DUl,l  ,l  ,l  ,DUl,l  ,l  ,DUl,l  ,l  ,l  ,DUl,l  ,l  ,
        /*C*/ u  ,uL ,uL ,UL ,DUL,DUL,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,
        /*G*/ u  ,uL ,DuL,uL ,DUL,DUl,DuL,L  ,l  ,Dul,l  ,Dul,l  ,Dul,l  ,Dul,l  ,Dul,l  ,l  ,Dul,l  ,Dul,l  ,Dul,l  ,Dul,l  ,
        /*A*/ u  ,DuL,uL ,DuL,DuL,ul ,Ul ,DUL,L  ,l  ,DUl,l  ,l  ,l  ,DUl,l  ,l  ,l  ,DUl,l  ,l  ,DUl,l  ,l  ,l  ,DUl,l  ,l  ,
        /*T*/ u  ,uL ,uL ,uL ,DuL,uL ,ul ,Ul ,DUL,DUL,DUl,DUl,DUl,Dul,Ul ,DUl,DUl,DUl,DUl,DUl,Dul,Ul ,DUl,DUl,DUl,DUl,DUl,DUl,
        /*C*/ u  ,uL ,uL ,uL ,DuL,Dul,uL ,ul ,DUl,DUL,Dul,l  ,Dul,l  ,ul ,DUl,Dul,L  ,Dul,Dul,l  ,l  ,l  ,Dul,l  ,Dul,l  ,Dul,
        /*G*/ u  ,uL ,DuL,uL ,DuL,Ul ,DuL,uL ,l  ,Dul,DUl,Dul,l  ,Dul,ul ,Dul,Ul ,Dul,DUl,DUl,Dul,Dul,Dul,Ul ,Dul,DUl,DUl,DUl,
        /*A*/ u  ,DuL,uL ,DuL,DuL,ul ,UL ,DuL,L  ,l  ,Dul,l  ,Dul,l  ,Dul,l  ,l  ,l  ,Dul,l  ,l  ,Dul,l  ,Dul,l  ,Dul,l  ,Dul,
        /*C*/ u  ,uL ,uL ,uL ,DuL,Dul,uL ,UL ,DuL,L  ,l  ,l  ,DUl,l  ,l  ,Dul,Dul,l  ,l  ,DUl,l  ,l  ,Dul,Dul,l  ,l  ,Dul,DUl,
        /*G*/ u  ,uL ,DuL,uL ,DuL,ul ,DuL,uL ,UL ,DUL,DUL,DUl,DUl,DUl,Dul,DUl,Dul,DUl,Dul,Ul ,DUl,L  ,DUl,l  ,l  ,l  ,DUl,l  ,
        /*A*/ u  ,DuL,uL ,DuL,DuL,ul ,uL ,DuL,uL ,DUL,Dul,L  ,Dul,l  ,Dul,l  ,Dul,l  ,Dul,l  ,Ul ,DUl,L  ,l  ,l  ,DUl,l  ,l  ,
        /*G*/ u  ,uL ,DuL,uL ,DuL,ul ,DuL,uL ,uL ,DuL,Ul ,DuL,L  ,Dul,l  ,Dul,l  ,l  ,l  ,l  ,Dul,Ul ,DUl,L  ,l  ,l  ,DUl,l  ,
        /*C*/ u  ,uL ,uL ,uL ,DuL,Dul,uL ,uL ,DuL,DuL,ul ,UL ,DUL,L  ,l  ,l  ,DUl,l  ,l  ,Dul,l  ,l  ,Ul ,DUl,L  ,l  ,l  ,DUl,
        /*G*/ u  ,uL ,DuL,uL ,DuL,ul ,DuL,uL ,uL ,DuL,ul ,DuL,UL ,DUL,L  ,Dul,l  ,DUl,l  ,l  ,Dul,l  ,Dul,Ul ,DUl,DUl,DUl,DUl,
        /*A*/ u  ,DuL,uL ,DuL,DuL,ul ,uL ,DuL,uL ,DuL,Dul,uL ,uL ,UL ,DUL,L  ,l  ,l  ,DUl,l  ,l  ,Dul,l  ,l  ,l  ,Dul,l  ,Dul,
        /*C*/ u  ,uL ,uL ,uL ,DuL,Dul,uL ,uL ,DuL,DuL,ul ,uL ,DuL,uL ,UL ,DUL,DUL,DUl,DUl,DUl,DUl,DUl,Dul,Dul,Dul,Ul ,Dul,DUl,
        /*T*/ u  ,uL ,uL ,uL ,DuL,uL ,ul ,uL ,uL ,DuL,uL ,ul ,uL ,uL ,uL ,DUL,DUl,DuL,L  ,l  ,Dul,l  ,l  ,l  ,Dul,l  ,l  ,DUl,
        /*A*/ u  ,DuL,uL ,DuL,DuL,Dul,uL ,DuL,uL ,DuL,Dul,uL ,ul ,uL ,DuL,DuL,ul ,Ul ,DUL,L  ,l  ,DUl,l  ,l  ,l  ,Dul,l  ,l  ,
        /*C*/ u  ,uL ,uL ,uL ,DuL,Dul,uL ,uL ,DuL,DuL,ul ,DuL,DuL,uL ,uL ,DuL,Dul,uL ,Ul ,DUL,L  ,l  ,l  ,DUl,l  ,l  ,DUl,DUl,
        /*G*/ u  ,uL ,DuL,uL ,DuL,ul ,DuL,uL ,uL ,DuL,ul ,DuL,uL ,DuL,uL ,DuL,uL ,ul ,ul ,UL ,DUL,L  ,DUl,l  ,l  ,l  ,DUl,l  ,
        /*T*/ u  ,uL ,uL ,uL ,DuL,uL ,ul ,uL ,uL ,DuL,uL ,ul ,uL ,uL ,uL ,DuL,Dul,DuL,ul ,uL ,UL ,DUL,DUL,DUl,DUl,DUl,DUl,DUl,
        /*A*/ u  ,DuL,uL ,DuL,DuL,Dul,uL ,DuL,uL ,DuL,Dul,uL ,ul ,uL ,DuL,DuL,ul ,ul ,DuL,uL ,uL ,DUL,DuL,l  ,Ul ,Dul,L  ,l  ,
        /*C*/ u  ,uL ,uL ,uL ,DuL,Dul,uL ,uL ,DuL,DuL,ul ,DuL,DuL,uL ,uL ,DuL,Dul,uL ,ul ,DuL,uL ,DuL,DUl,DuL,L  ,Ul ,DUl,DUl,
        /*G*/ u  ,uL ,DuL,uL ,DuL,ul ,DuL,uL ,uL ,DuL,ul ,DuL,uL ,DuL,uL ,DuL,uL ,ul ,ul ,uL ,DuL,DuL,Dul,DUL,Dul,uL ,DUl,DUL,
        /*A*/ u  ,DuL,uL ,DuL,DuL,ul ,uL ,DuL,uL ,DuL,Dul,uL ,uL ,uL ,DuL,DuL,Dul,uL ,Dul,uL ,uL ,DuL,L  ,Dul,DUl,Dul,L  ,Dul
        }
    };
}();

static auto dna4_match_4_mismatch_5_gap_1_open_10_part_04 = []()
{
    using seqan3::operator""_dna4;
    return alignment_fixture
    {
        "AGATCGACTAGCG"_dna4,
        "G"_dna4,
        align_config_dna_score,
        -18,
        "AGATCGACTAGCG",
        "------------G",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 13u,
        /*.sequence2_end_position = */ 1u,
        std::vector
        {
        //      e,  A,  G,  A,  T,  C,  G,  A,  C,  T,  A,  G,  C,  G,
        /*e*/ 0  ,-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,-22,-23,
        /*G*/ -11,-5 ,-7 ,-17,-18,-19,-11,-21,-22,-23,-24,-16,-26,-18
        },
        std::vector
        {
        //      e,  A,  G,  A,  T,  C,  G,  A,  C,  T,  A,  G,  C,  G,
        /*e*/ N  ,L  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*G*/ U  ,DUL,DUL,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl
        }
    };
}();

static auto dna4_match_4_mismatch_5_gap_1_open_10_part_05 = []()
{
    using seqan3::operator""_dna4;
    return alignment_fixture
    {
        "A"_dna4,
        "AGATCGACTAGCG"_dna4,
        align_config_dna_score,
        -18,
        "A------------",
        "AGATCGACTAGCG",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 1u,
        /*.sequence2_end_position = */ 13u,
        std::vector
        {
        //       e,  A,
        /*e*/  0  ,-11,
        /*A*/  -11,4  ,
        /*G*/  -12,-7 ,
        /*A*/  -13,-8 ,
        /*T*/  -14,-9 ,
        /*C*/  -15,-10,
        /*G*/  -16,-11,
        /*A*/  -17,-12,
        /*C*/  -18,-13,
        /*T*/  -19,-14,
        /*A*/  -20,-15,
        /*G*/  -21,-16,
        /*C*/  -22,-17,
        /*G*/  -23,-18
        },
        std::vector
        {
        //       e,  A,
        /*e*/  N  ,L  ,
        /*A*/  U  ,DUL,
        /*G*/  u  ,UL ,
        /*A*/  u  ,DuL,
        /*T*/  u  ,uL ,
        /*C*/  u  ,uL ,
        /*G*/  u  ,uL ,
        /*A*/  u  ,DuL,
        /*C*/  u  ,uL ,
        /*T*/  u  ,uL ,
        /*A*/  u  ,DuL,
        /*G*/  u  ,uL ,
        /*C*/  u  ,uL ,
        /*G*/  u  ,uL
        }
    };
}();

static auto dna4_match_4_mismatch_5_gap_1_open_10_seq1_empty = []()
{
    using seqan3::operator""_dna4;
    return alignment_fixture
    {
        ""_dna4,
        "AGATCGACT"_dna4,
        align_config_dna_score,
        -19,
        "---------",
        "AGATCGACT",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 0u,
        /*.sequence2_end_position = */ 9u,
        std::vector
        {
        //       e,
        /*e*/  0  ,
        /*A*/  -11,
        /*G*/  -12,
        /*A*/  -13,
        /*T*/  -14,
        /*C*/  -15,
        /*G*/  -16,
        /*A*/  -17,
        /*C*/  -18,
        /*T*/  -19
        },
        std::vector
        {
        //     e,
        /*e*/  N,
        /*A*/  U,
        /*G*/  u,
        /*A*/  u,
        /*T*/  u,
        /*C*/  u,
        /*G*/  u,
        /*A*/  u,
        /*C*/  u,
        /*T*/  u
        }
    };
}();

static auto dna4_match_4_mismatch_5_gap_1_open_10_seq2_empty = []()
{
    using seqan3::operator""_dna4;
    return alignment_fixture
    {
        "AGATCGACT"_dna4,
        ""_dna4,
        align_config_dna_score,
        -19,
        "AGATCGACT",
        "---------",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 9u,
        /*.sequence2_end_position = */ 0u,
        std::vector
        {
        //      e,  A,  G,  A,  T,  C,  G,  A,  C,  T,
        /*e*/ 0  ,-11,-12,-13,-14,-15,-16,-17,-18,-19,
        },
        std::vector
        {
        //    e,A,G,A,T,C,G,A,C,T
        /*e*/ N,L,l,l,l,l,l,l,l,l
        }
    };
}();

static auto dna4_match_4_mismatch_5_gap_1_open_10_both_empty = []()
{
    using seqan3::operator""_dna4;
    return alignment_fixture
    {
        ""_dna4,
        ""_dna4,
        align_config_dna_score,
        0,
        "",
        "",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 0u,
        /*.sequence2_end_position = */ 0u,
        std::vector
        {
        //    e
        /*e*/ 0
        },
        std::vector
        {
        //    e
        /*e*/ N
        }
    };
}();

// [ISSUE #3043] Wrong alignment score https://github.com/seqan/seqan3/issues/3043
static auto issue_3043 = []()
{
    using seqan3::operator""_dna4;

    auto make_scheme = [] ()
    {
        seqan3::nucleotide_scoring_scheme<int32_t> scheme{};
        scheme.score('A'_dna4, 'A'_dna4) = 145;
        scheme.score('A'_dna4, 'T'_dna4) = -144;
        scheme.score('A'_dna4, 'G'_dna4) = -153;
        scheme.score('A'_dna4, 'C'_dna4) = -152;

        scheme.score('T'_dna4, 'A'_dna4) = -144;
        scheme.score('T'_dna4, 'T'_dna4) = 144;
        scheme.score('T'_dna4, 'G'_dna4) = -143;
        scheme.score('T'_dna4, 'C'_dna4) = -140;

        scheme.score('G'_dna4, 'A'_dna4) = -153;
        scheme.score('G'_dna4, 'T'_dna4) = -143;
        scheme.score('G'_dna4, 'G'_dna4) = 144;
        scheme.score('G'_dna4, 'C'_dna4) = -145;

        scheme.score('C'_dna4, 'A'_dna4) = -152;
        scheme.score('C'_dna4, 'T'_dna4) = -140;
        scheme.score('C'_dna4, 'G'_dna4) = -145;
        scheme.score('C'_dna4, 'C'_dna4) = 144;

        return scheme;
    };

    return alignment_fixture
    {
        "GGCAAGAA"_dna4,
        "CGAAGC"_dna4,
        seqan3::align_cfg::method_global{} | seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{-52},
                                                                                seqan3::align_cfg::extension_score{-58}}
                                           | seqan3::align_cfg::scoring_scheme{make_scheme()},
        74,
        "GGCAAGAA--",
        "--C--GAAGC",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 8u,
        /*.sequence2_end_position = */ 6u,
        std::vector
        {
        //    e
        /*e*/    0,-110,-168,-226,-284,-342,-400,-458,-516,
              -110,-145,-255, -24,-134,-192,-250,-308,-366,
              -168,  34,  -1,-111,-169,-227, -48,-158,-216,
              -226, -76,-111,-153,  34, -24,-134,  97, -13,
              -284,-134,-169,-250,  -8, 179,  69,  11, 242,
              -342,-140,  10,-100,-118,  69, 323, 213, 155,
              -400,-250,-100, 154,  44,  11, 213, 171,  74
        },
        std::vector
        {
        //    e,G  ,G  ,C  ,A  ,A  ,G  ,A  ,A  ,
        /*e*/ N,L  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*C*/ U,DUL,DUL,DUl,L  ,l  ,l  ,l  ,l  ,
        /*G*/ u,DUL,DuL,L  ,l  ,l  ,DUl,L  ,l  ,
        /*A*/ u,UL ,UL ,DuL,DUL,DUL,l  ,DUl,DUL,
        /*A*/ u,uL ,uL ,uL ,DUl,DUL,L  ,DUl,DUl,
        /*G*/ u,DuL,DuL,L  ,Ul ,Ul ,DUL,L  ,l  ,
        /*C*/ u,uL ,UL ,DUL,L  ,ul ,Ul ,DUL,uL
        }
    };
}();

// ----------------------------------------------------------------------------
// alignment fixtures using amino acid alphabet
// ----------------------------------------------------------------------------

inline constexpr auto config_blosum62_scheme =
    seqan3::align_cfg::scoring_scheme{seqan3::aminoacid_scoring_scheme{seqan3::aminoacid_similarity_matrix::blosum62}};

static auto aa27_blosum62_gap_1_open_10 = [] ()
{
    using seqan3::operator""_aa27;
    return alignment_fixture
    {
        "FNQSAEYPDISHCGVMQLKWRATLGT"_aa27,
        "EIKSDVLLHRWSMKNPGNILMIDVGMQVAESYFAT"_aa27,
        align_config | config_blosum62_scheme,
        -26,
        "--------FNQSAEYP-DISHCGVMQLKWRATLGT",
        "EIKSDVLLHRWSMKNPGNILMIDVGMQVAESYFAT",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 26u,
        /*.sequence2_end_position = */ 35u,
        std::vector
        {
        //    e  ,F  ,N  ,Q  ,S  ,A  ,E  ,Y  ,P  ,D  ,I  ,S  ,H  ,C  ,G  ,V  ,M  ,Q  ,L  ,K  ,W  ,R  ,A  ,T  ,L  ,G  ,T  ,
        /*e*/ 0  ,-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,-22,-23,-24,-25,-26,-27,-28,-29,-30,-31,-32,-33,-34,-35,-36,
        /*E*/ -11,-3 ,-11,-10,-13,-15,-10,-18,-18,-16,-22,-20,-21,-25,-25,-26,-27,-24,-30,-27,-32,-30,-32,-33,-36,-36,-36,
        /*I*/ -12,-11,-6 ,-14,-12,-14,-18,-11,-21,-21,-12,-23,-23,-22,-26,-22,-25,-29,-22,-31,-30,-33,-31,-33,-31,-37,-37,
        /*K*/ -13,-15,-11,-5 ,-14,-13,-13,-19,-12,-21,-22,-12,-23,-24,-24,-26,-23,-24,-29,-17,-28,-28,-30,-31,-32,-33,-34,
        /*S*/ -14,-15,-14,-11,-1 ,-12,-13,-14,-15,-12,-17,-18,-13,-20,-21,-22,-23,-23,-25,-26,-20,-28,-27,-29,-31,-32,-32,
        /*D*/ -15,-17,-14,-14,-11,-3 ,-10,-15,-15,-9 ,-15,-17,-19,-16,-21,-23,-24,-23,-26,-26,-28,-22,-30,-28,-32,-32,-33,
        /*V*/ -16,-16,-20,-16,-13,-11,-5 ,-11,-17,-18,-6 ,-17,-18,-19,-19,-17,-22,-23,-22,-25,-26,-27,-22,-29,-27,-31,-32,
        /*L*/ -17,-16,-19,-19,-14,-14,-14,-6 ,-14,-18,-16,-8 ,-19,-19,-21,-18,-15,-24,-19,-24,-27,-28,-28,-23,-25,-31,-32,
        /*L*/ -18,-17,-19,-20,-15,-15,-17,-15,-9 ,-18,-16,-18,-11,-20,-23,-20,-16,-17,-20,-21,-26,-29,-29,-29,-19,-29,-31,
        /*H*/ -19,-19,-16,-19,-16,-17,-15,-15,-17,-10,-19,-17,-10,-14,-22,-23,-22,-16,-20,-21,-23,-26,-30,-31,-30,-21,-31,
        /*R*/ -20,-22,-19,-15,-17,-17,-17,-17,-17,-19,-13,-20,-17,-13,-16,-25,-24,-21,-18,-18,-24,-18,-27,-30,-31,-32,-22,
        /*W*/ -21,-19,-25,-21,-18,-19,-20,-15,-21,-21,-21,-16,-22,-19,-15,-19,-26,-26,-23,-21,-7 ,-18,-19,-20,-21,-22,-23,
        /*S*/ -22,-23,-18,-24,-17,-17,-19,-21,-16,-21,-22,-17,-17,-23,-19,-17,-20,-26,-28,-23,-18,-8 ,-17,-18,-21,-21,-21,
        /*M*/ -23,-22,-25,-18,-20,-18,-19,-20,-23,-19,-20,-23,-19,-18,-26,-18,-12,-20,-24,-25,-19,-19,-9 ,-18,-16,-22,-22,
        /*K*/ -24,-26,-22,-24,-18,-21,-17,-21,-21,-24,-22,-20,-24,-22,-20,-28,-19,-11,-22,-19,-20,-17,-20,-10,-20,-18,-23,
        /*N*/ -25,-27,-20,-22,-22,-20,-21,-19,-23,-20,-25,-21,-19,-27,-22,-23,-24,-19,-14,-22,-21,-20,-19,-20,-13,-20,-18,
        /*P*/ -26,-28,-29,-21,-23,-23,-21,-24,-12,-23,-23,-25,-23,-22,-28,-24,-25,-23,-22,-15,-22,-22,-21,-20,-23,-15,-21,
        /*G*/ -27,-29,-28,-29,-21,-23,-25,-24,-23,-13,-24,-23,-26,-26,-16,-27,-26,-24,-26,-24,-17,-23,-22,-23,-24,-17,-17,
        /*N*/ -28,-30,-23,-28,-25,-23,-23,-27,-24,-22,-16,-23,-22,-29,-26,-19,-27,-25,-27,-26,-24,-17,-24,-22,-26,-24,-17,
        /*I*/ -29,-28,-33,-26,-26,-26,-26,-24,-25,-25,-18,-18,-26,-23,-28,-23,-18,-26,-23,-28,-25,-25,-18,-25,-20,-28,-25,
        /*L*/ -30,-29,-31,-32,-27,-27,-29,-27,-26,-26,-23,-20,-21,-27,-27,-27,-21,-20,-22,-25,-26,-26,-26,-19,-21,-24,-29,
        /*M*/ -31,-30,-31,-31,-28,-28,-29,-30,-27,-27,-25,-24,-22,-22,-30,-26,-22,-21,-18,-23,-26,-27,-27,-27,-17,-24,-25,
        /*I*/ -32,-31,-33,-34,-29,-29,-31,-30,-28,-28,-23,-27,-27,-23,-26,-27,-25,-25,-19,-21,-26,-28,-28,-28,-25,-21,-25,
        /*D*/ -33,-35,-30,-33,-30,-31,-27,-32,-29,-22,-31,-23,-28,-30,-24,-29,-30,-25,-29,-20,-25,-28,-29,-29,-29,-26,-22,
        /*V*/ -34,-34,-38,-32,-31,-30,-33,-28,-30,-30,-19,-30,-26,-29,-33,-20,-28,-31,-24,-31,-23,-28,-28,-29,-28,-32,-26,
        /*G*/ -35,-37,-34,-37,-32,-31,-32,-34,-30,-31,-30,-19,-30,-29,-23,-31,-23,-30,-32,-26,-31,-25,-28,-30,-31,-22,-33,
        /*M*/ -36,-35,-39,-34,-33,-33,-33,-33,-32,-32,-30,-30,-21,-31,-32,-22,-26,-23,-28,-33,-27,-32,-26,-29,-28,-33,-23,
        /*Q*/ -37,-39,-35,-34,-34,-34,-31,-34,-33,-32,-32,-30,-30,-24,-33,-33,-22,-21,-25,-27,-33,-26,-33,-27,-31,-30,-34,
        /*V*/ -38,-38,-42,-37,-35,-34,-36,-32,-34,-34,-29,-32,-33,-31,-27,-29,-32,-24,-20,-27,-30,-33,-26,-33,-26,-34,-30,
        /*A*/ -39,-40,-40,-41,-36,-31,-35,-38,-33,-35,-34,-28,-34,-33,-31,-27,-30,-33,-25,-21,-30,-31,-29,-26,-34,-26,-34,
        /*E*/ -40,-42,-40,-38,-37,-37,-26,-37,-36,-31,-35,-34,-28,-37,-35,-33,-29,-28,-32,-24,-24,-30,-32,-30,-29,-36,-27,
        /*S*/ -41,-42,-41,-40,-34,-36,-37,-28,-37,-36,-33,-31,-35,-29,-37,-37,-34,-29,-30,-32,-27,-25,-29,-31,-32,-29,-35,
        /*Y*/ -42,-38,-44,-42,-39,-36,-38,-30,-31,-38,-37,-35,-29,-37,-32,-38,-37,-35,-30,-32,-30,-29,-27,-31,-32,-35,-31,
        /*F*/ -43,-36,-41,-45,-40,-41,-39,-35,-34,-34,-38,-37,-36,-31,-40,-33,-38,-37,-35,-33,-31,-33,-31,-29,-31,-35,-37,
        /*A*/ -44,-45,-38,-42,-41,-36,-40,-41,-36,-36,-35,-37,-39,-36,-31,-40,-34,-38,-36,-36,-36,-32,-29,-31,-30,-31,-35,
        /*T*/ -45,-46,-45,-39,-41,-41,-37,-42,-41,-37,-37,-34,-39,-40,-38,-31,-40,-35,-37,-37,-38,-37,-32,-24,-32,-32,-26
        },
        std::vector
        {
        //    e  ,F  ,N  ,Q  ,S  ,A  ,E  ,Y  ,P  ,D  ,I  ,S  ,H  ,C  ,G  ,V  ,M  ,Q  ,L  ,K  ,W  ,R  ,A  ,T  ,L  ,G  ,T  ,
        /*e*/ N  ,L  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*E*/ U  ,DUL,DUL,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,l  ,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,
        /*I*/ u  ,DUL,DUL,DUL,DUl,DUl,DUl,Dul,DUl,DUl,Dul,L  ,DUl,Dul,l  ,Dul,Dul,l  ,Dul,l  ,Dul,l  ,DUl,DUl,Dul,l  ,DUl,
        /*K*/ u  ,DuL,DUl,DuL,DUL,DUl,Dul,l  ,Dul,l  ,l  ,Dul,L  ,l  ,Dul,l  ,DUl,Dul,l  ,Dul,L  ,Dul,l  ,l  ,l  ,Dul,l  ,
        /*S*/ u  ,DuL,Dul,DUL,DuL,L  ,Dul,l  ,l  ,Dul,l  ,DUl,Dul,l  ,l  ,l  ,l  ,DUl,l  ,l  ,DUl,l  ,DUl,DUl,l  ,DUl,DUl,
        /*D*/ u  ,DuL,Dul,DuL,DUL,DUL,DuL,l  ,Dul,DUl,Dul,Dul,DUl,DUl,DUl,l  ,l  ,DUl,l  ,Dul,l  ,DUl,DUl,DUl,l  ,DUl,DUl,
        /*V*/ u  ,DuL,DuL,Dul,uL ,DUL,DUL,DuL,Dul,DUl,Dul,DuL,l  ,l  ,DUl,Dul,Dul,l  ,Dul,l  ,l  ,l  ,Dul,l  ,Dul,l  ,Dul,
        /*L*/ u  ,DuL,DuL,ul ,ul ,DuL,DUL,DUL,DuL,l  ,DUl,Dul,L  ,Dul,l  ,DUl,DUl,DUl,DUl,Dul,Dul,Dul,DUl,Dul,DUl,DUl,DUl,
        /*L*/ u  ,DuL,DuL,ul ,ul ,DuL,DuL,DUl,DUL,DuL,Dul,DUl,Dul,DuL,Dul,Dul,DUl,Dul,DUl,Dul,Dul,Dul,Dul,DUl,DUl,DUL,l  ,
        /*H*/ u  ,DuL,DuL,DuL,ul ,DuL,Dul,DuL,DUL,Dul,uL ,Dul,DUl,DuL,Dul,l  ,Dul,DUl,Dul,DUl,Dul,Dul,l  ,Dul,Ul ,DUl,DUL,
        /*R*/ u  ,DuL,Dul,DuL,uL ,Dul,Dul,DuL,DuL,DUL,Dul,DuL,DUl,DUl,DuL,Dul,Dul,DUl,DUl,DUl,DUL,Dul,DuL,l  ,ul ,DUl,DUl,
        /*W*/ u  ,DuL,uL ,Dul,Dul,uL ,Dul,Dul,DuL,Dul,ul ,Dul,DuL,DUl,DUl,DuL,Dul,Dul,DUl,DUl,Dul,L  ,l  ,l  ,l  ,l  ,l  ,
        /*S*/ u  ,DuL,Dul,uL ,Dul,DuL,DuL,ul ,Dul,DuL,ul ,Dul,DuL,DuL,DUl,DUl,DuL,Dul,Dul,Dul,Ul ,DUL,DUL,DUl,l  ,DUl,Dul,
        /*M*/ u  ,DuL,DuL,Dul,uL ,Dul,DuL,Dul,Dul,Dul,DuL,Dul,Dul,DuL,DuL,DUl,DuL,DuL,Dul,l  ,ul ,DUl,DUl,DUL,DUl,l  ,DUl,
        /*K*/ u  ,DuL,Dul,DuL,Dul,DuL,Dul,DuL,Dul,Dul,Dul,Dul,DuL,Dul,Dul,DuL,DUl,DuL,DuL,Dul,ul ,Dul,DUl,DUl,DUL,Dul,Dul,
        /*N*/ u  ,DuL,Dul,DuL,ul ,Dul,DuL,Dul,DuL,Dul,uL ,Dul,DuL,DuL,Dul,Dul,ul ,DUl,DuL,DUL,ul ,Dul,Dul,DUl,Dul,DUL,Dul,
        /*P*/ u  ,uL ,Dul,Dul,DuL,Dul,Dul,DuL,Dul,L  ,Dul,l  ,Dul,Dul,l  ,Dul,Dul,ul ,DUl,Dul,uL ,ul ,Dul,Dul,DUl,Dul,DUL,
        /*G*/ u  ,DuL,Dul,uL ,Dul,DuL,Dul,Dul,Ul ,DuL,L  ,Dul,l  ,Dul,Dul,L  ,ul ,ul ,ul ,DUl,Dul,uL ,Dul,Dul,Dul,DUl,DuL,
        /*N*/ u  ,DuL,Dul,DuL,ul ,Dul,DuL,DuL,ul ,DUL,DuL,DuL,Dul,Dul,DUl,Dul,uL ,ul ,Dul,Dul,ul ,Dul,uL ,Dul,Dul,Dul,DUl,
        /*I*/ u  ,DuL,DuL,Dul,uL ,DuL,DuL,DuL,uL ,ul ,DUL,DuL,DuL,Dul,ul ,DUl,Dul,uL ,Dul,ul ,ul ,ul ,Dul,DuL,Dul,ul ,DUl,
        /*L*/ u  ,DuL,DuL,ul ,ul ,DuL,DuL,Dul,uL ,uL ,DuL,DUL,DuL,Dul,Dul,Dul,Dul,DuL,DuL,Dul,ul ,ul ,Dul,Dul,DuL,Dul,Dul,
        /*M*/ u  ,DuL,DuL,Dul,uL ,DuL,DuL,Dul,ul ,uL ,DuL,DuL,DuL,DuL,DuL,Dul,Dul,DuL,DuL,DuL,Dul,Dul,Dul,Dul,Dul,DuL,Dul,
        /*I*/ u  ,DuL,DuL,Dul,ul ,DuL,DuL,Dul,uL ,uL ,DuL,DuL,Dul,DUl,DuL,Dul,Dul,DuL,DUL,DuL,Dul,ul ,Dul,Dul,DUl,Dul,DuL,
        /*D*/ u  ,DuL,Dul,DuL,ul ,DuL,Dul,uL ,ul ,Dul,DuL,Dul,DuL,Dul,Dul,DuL,Dul,Dul,DuL,Dul,DuL,Dul,ul ,Dul,ul ,Dul,Dul,
        /*V*/ u  ,DuL,DuL,Dul,uL ,DuL,DuL,Dul,uL ,ul ,Dul,L  ,Dul,Dul,Dul,Dul,DuL,ul ,Dul,DUl,Dul,DuL,Dul,Dul,Dul,Dul,Dul,
        /*G*/ u  ,DuL,Dul,uL ,Dul,DuL,DuL,ul ,Dul,DuL,Ul ,DuL,L  ,Dul,Dul,Ul ,Dul,DuL,ul ,Dul,ul ,Dul,DuL,Dul,ul ,Dul,L  ,
        /*M*/ u  ,DuL,DuL,Dul,uL ,DuL,DuL,DuL,uL ,uL ,DuL,UL ,DuL,DuL,DUl,Dul,DUL,Dul,DuL,Dul,Dul,Dul,Dul,DuL,Dul,Ul ,Dul,
        /*Q*/ u  ,DuL,Dul,DuL,DuL,DuL,DuL,DuL,ul ,Dul,uL ,DuL,DUL,DuL,DuL,ul ,Dul,DuL,DuL,Dul,ul ,Dul,Dul,Dul,Dul,Dul,DUl,
        /*V*/ u  ,DuL,DuL,Dul,uL ,DuL,DuL,Dul,uL ,ul ,Dul,uL ,Dul,DUl,DuL,DuL,DUl,DUl,DuL,DuL,Dul,l  ,Dul,Dul,Dul,Dul,Dul,
        /*A*/ u  ,DuL,Dul,uL ,Dul,DuL,DuL,Dul,Dul,uL ,ul ,DuL,DuL,Dul,Dul,Dul,DuL,Dul,DUl,DuL,DuL,Dul,Dul,Dul,Dul,Dul,DuL,
        /*E*/ u  ,DuL,Dul,DuL,uL ,DuL,DuL,DuL,ul ,Dul,ul ,Dul,Dul,uL ,Dul,Dul,Dul,DuL,uL ,DUl,DuL,DuL,Dul,Dul,Dul,Dul,Dul,
        /*S*/ u  ,DuL,Dul,DuL,DuL,DuL,DUl,Dul,uL ,Dul,Dul,Dul,DuL,Dul,DuL,Dul,Dul,Dul,DuL,Dul,DUl,DuL,DuL,Dul,Dul,Dul,Dul,
        /*Y*/ u  ,DuL,DuL,Dul,ul ,DuL,DuL,DUl,DuL,ul ,Dul,Dul,Dul,DuL,Dul,Dul,ul ,Dul,Dul,DuL,Dul,DUL,DuL,DuL,Dul,Dul,Dul,
        /*F*/ u  ,DuL,DuL,ul ,ul ,Dul,Dul,DuL,DuL,DuL,DuL,ul ,Dul,DuL,DuL,Dul,Dul,ul ,Dul,DuL,DuL,DuL,DUl,DuL,DuL,Dul,Dul,
        /*A*/ u  ,DuL,Dul,DuL,ul ,Dul,uL ,Dul,Dul,DuL,DuL,DuL,Dul,Dul,DuL,DuL,Dul,ul ,ul ,Dul,Dul,DuL,DuL,DuL,Dul,DuL,Dul,
        /*T*/ u  ,DuL,Dul,DuL,DuL,Dul,Dul,DuL,ul ,Dul,DuL,DuL,DuL,Dul,DUl,Dul,uL ,Dul,ul ,Dul,Dul,Dul,Dul,DuL,DuL,Dul,Dul
        }
    };
}();

static auto aa27_blosum62_gap_1_open_10_small = [] ()
{
    using seqan3::operator""_aa27;
    return alignment_fixture
    {
        "RKFCYMD"_aa27,
        "GAYQW"_aa27,
        align_config | config_blosum62_scheme,
        -11,
        "RKFCYMD",
        "G--AYQW",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 7u,
        /*.sequence2_end_position = */ 5u,
        std::vector
        {
        //    e  ,R  ,K  ,F  ,C  ,Y  ,M  ,D  ,
        /*e*/  0 ,-11,-12,-13,-14,-15,-16,-17,
        /*G*/ -11,-2 ,-13,-14,-15,-16,-17,-17,
        /*A*/ -12,-12,-3 ,-14,-14,-16,-17,-18,
        /*Y*/ -13,-14,-14,0  ,-11,-7 ,-13,-14,
        /*Q*/ -14,-12,-13,-11,-3 ,-12,-7 ,-13,
        /*W*/ -15,-16,-15,-12,-13,-1 ,-12,-11
        },
        std::vector
        {
        //    e  ,R  ,K  ,F  ,C  ,Y  ,M  ,D  ,
        /*e*/ N  ,L  ,l  ,l  ,l  ,l  ,l  ,l  ,
        /*G*/ U  ,DUL,DUL,l  ,l  ,l  ,l  ,DUl,
        /*A*/ u  ,DUL,DuL,L  ,Dul,l  ,Dul,l  ,
        /*Y*/ u  ,DuL,DUl,DUL,L  ,DUl,l  ,l  ,
        /*Q*/ u  ,DuL,DuL,Ul ,DUL,DUL,DUl,DUl,
        /*W*/ u  ,uL ,Dul,DuL,DUL,Dul,L  ,DUl
        }
    };
}();

static auto aa27_blosum62_gap_1_open_10_empty_first = [] ()
{
    using seqan3::operator""_aa27;
    return alignment_fixture
    {
        "PPAMDYIRPW"_aa27,
        ""_aa27,
        align_config | config_blosum62_scheme,
        -20,
        "PPAMDYIRPW",
        "----------",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 10u,
        /*.sequence2_end_position = */ 0u,
        std::vector
        {
        //    e  ,P  ,P  ,A  ,M  ,D  ,Y  ,I  ,R  ,P  ,W  ,
        /*e*/ 0  ,-11,-12,-13,-14,-15,-16,-17,-18,-19,-20
        },
        std::vector
        {
        //    e,P,P,A,M,D,Y,I,R,P,W,
        /*e*/ N,L,l,l,l,l,l,l,l,l,l
        }
    };
}();

static auto aa27_blosum62_gap_1_open_10_empty_second = [] ()
{
    using seqan3::operator""_aa27;
    return alignment_fixture
    {
        ""_aa27,
        "PPAMDYIRPW"_aa27,
        align_config | config_blosum62_scheme,
        -20,
        "----------",
        "PPAMDYIRPW",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 0u,
        /*.sequence2_end_position = */ 10u,
        std::vector
        {
        //    e
        /*e*/ 0  ,
        /*P*/ -11,
        /*P*/ -12,
        /*A*/ -13,
        /*M*/ -14,
        /*D*/ -15,
        /*Y*/ -16,
        /*I*/ -17,
        /*R*/ -18,
        /*P*/ -19,
        /*W*/ -20
        },
        std::vector
        {
        //    e
        /*e*/ N,
        /*P*/ U,
        /*P*/ u,
        /*A*/ u,
        /*M*/ u,
        /*D*/ u,
        /*Y*/ u,
        /*I*/ u,
        /*R*/ u,
        /*P*/ u,
        /*W*/ u
        }
    };
}();

static auto aa27_blosum62_gap_1_open_10_empty_both = [] ()
{
    using seqan3::operator""_aa27;
    return alignment_fixture
    {
        ""_aa27,
        ""_aa27,
        align_config | config_blosum62_scheme,
        0,
        "",
        "",
        /*.sequence1_begin_position = */ 0u,
        /*.sequence2_begin_position = */ 0u,
        /*.sequence1_end_position = */ 0u,
        /*.sequence2_end_position = */ 0u,
        std::vector
        {
        //     e
        /*e*/  0
        },
        std::vector
        {
        //    e
        /*e*/ N
        }
    };
}();
// clang-format on

} // namespace seqan3::test::alignment::fixture::global::affine::unbanded
