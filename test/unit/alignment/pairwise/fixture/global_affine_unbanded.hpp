// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
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

using seqan3::operator""_dna4;

namespace seqan3::test::alignment::fixture::global::affine::unbanded
{

inline constexpr auto align_config = seqan3::align_cfg::mode{seqan3::global_alignment}
                                   | seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1},
                                                            seqan3::gap_open_score{-10}}};
inline constexpr auto align_config_dna_score = align_config
                                             | seqan3::align_cfg::scoring{
                                                    seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                      seqan3::mismatch_score{-5}}};

// Naming is used to disambiguate configs for global alignment with different scoring parameters.
// `part_XX` indicates that the config is the same but the sequences differ in their lengths.
// This is helpful to collect different test fixtures to test vectorised alignments where we can
// only apply alignments with the same config.
static auto dna4_match_4_mismatch_5_gap_1_open_10_part_01 = []()
{
    using seqan3::operator""_dna4;
    return seqan3::test::alignment::fixture::alignment_fixture
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
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}},
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{16u}, seqan3::detail::row_index_type{9u}},
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
    return seqan3::test::alignment::fixture::alignment_fixture
    {
        "ACGTACGTA"_dna4,
        "AACCGGTTAACCGGTT"_dna4,
        align_config_dna_score,
        -18,
        "ACGTAC----------GTA",
        "A---ACCGGTTAACCGGTT",
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}},
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{9u}, seqan3::detail::row_index_type{16u}},
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
    return seqan3::test::alignment::fixture::alignment_fixture
    {
        "AGATCGACTAGCGAGCTACGAGCTAGC"_dna4,
        "AGACGATCGACGAGCGACTACGTACGA"_dna4,
        align_config_dna_score,
        26,
        "A---GATCGACTAGCGAGCTACGAGCTA-GC",
        "AGACGATCGACGAGCGA-CTACG---TACGA",
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}},
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{27u}, seqan3::detail::row_index_type{27u}},
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
    return seqan3::test::alignment::fixture::alignment_fixture
    {
        "AGATCGACTAGCG"_dna4,
        "G"_dna4,
        align_config_dna_score,
        -18,
        "AGATCGACTAGCG",
        "------------G",
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}},
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{13u}, seqan3::detail::row_index_type{1u}},
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
    return seqan3::test::alignment::fixture::alignment_fixture
    {
        "A"_dna4,
        "AGATCGACTAGCG"_dna4,
        align_config_dna_score,
        -18,
        "A------------",
        "AGATCGACTAGCG",
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}},
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{1u}, seqan3::detail::row_index_type{13u}},
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
    return seqan3::test::alignment::fixture::alignment_fixture
    {
        ""_dna4,
        "AGATCGACT"_dna4,
        align_config_dna_score,
        -19,
        "---------",
        "AGATCGACT",
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}},
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{9u}},
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
    return seqan3::test::alignment::fixture::alignment_fixture
    {
        "AGATCGACT"_dna4,
        ""_dna4,
        align_config_dna_score,
        -19,
        "AGATCGACT",
        "---------",
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}},
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{9u}, seqan3::detail::row_index_type{0u}},
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
    return seqan3::test::alignment::fixture::alignment_fixture
    {
        ""_dna4,
        ""_dna4,
        align_config_dna_score,
        0,
        "",
        "",
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}},
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}},
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

} // namespace seqan3::test::alignment::fixture::global::affine::unbanded
