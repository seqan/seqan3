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
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>

#include "alignment_fixture.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_rna5;

namespace seqan3::test::alignment::fixture::local::affine::unbanded
{

inline constexpr auto align_config = seqan3::align_cfg::mode{seqan3::local_alignment} |
                                     seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1},
                                                                               seqan3::gap_open_score{-10}}};

// Local alignment with mismatch.
static auto dna4_01 = []()
{
    using seqan3::detail::column_index_type;
    using seqan3::detail::row_index_type;

    return alignment_fixture
    {
        // score: 11 (4 matches, 1 mismatch)
        // alignment:
        // GTTTA
        // || ||
        // GTCTA
        "AACCGGTTTAACCGGTT"_dna4,
        "ACGTCTACGTA"_dna4,
        align_config | seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                    seqan3::mismatch_score{-5}}},
        11,
        "GTTTA",
        "GTCTA",
        seqan3::alignment_coordinate{column_index_type{5u}, row_index_type{2u}},
        seqan3::alignment_coordinate{column_index_type{10u}, row_index_type{7u}},
        std::vector
        {
        //     e, A, A, C, C, G, G, T, T, T, A, A, C, C, G, G, T, T
        /*e*/ 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
        /*A*/ 0 ,4 ,4 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4 ,4 ,0 ,0 ,0 ,0 ,0 ,0 ,
        /*C*/ 0 ,0 ,0 ,8 ,4 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,8 ,4 ,0 ,0 ,0 ,0 ,
        /*G*/ 0 ,0 ,0 ,0 ,3 ,8 ,4 ,0 ,0 ,0 ,0 ,0 ,0 ,3 ,8 ,4 ,0 ,0 ,
        /*T*/ 0 ,0 ,0 ,0 ,0 ,0 ,3 ,8 ,4 ,4 ,0 ,0 ,0 ,0 ,0 ,3 ,8 ,4 ,
        /*C*/ 0 ,0 ,0 ,4 ,4 ,0 ,0 ,0 ,3 ,0 ,0 ,0 ,4 ,4 ,0 ,0 ,0 ,3 ,
        /*T*/ 0 ,0 ,0 ,0 ,0 ,0 ,0 ,4 ,4 ,7 ,0 ,0 ,0 ,0 ,0 ,0 ,4 ,4 ,
        /*A*/ 0 ,4 ,4 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,11,4 ,0 ,0 ,0 ,0 ,0 ,0 ,
        /*C*/ 0 ,0 ,0 ,8 ,4 ,0 ,0 ,0 ,0 ,0 ,0 ,6 ,8 ,4 ,0 ,0 ,0 ,0 ,
        /*G*/ 0 ,0 ,0 ,0 ,3 ,8 ,4 ,0 ,0 ,0 ,0 ,0 ,1 ,3 ,8 ,4 ,0 ,0 ,
        /*T*/ 0 ,0 ,0 ,0 ,0 ,0 ,3 ,8 ,4 ,4 ,0 ,0 ,0 ,0 ,0 ,3 ,8 ,4 ,
        /*A*/ 0 ,4 ,4 ,0 ,0 ,0 ,0 ,0 ,3 ,0 ,8 ,4 ,0 ,0 ,0 ,0 ,0 ,3
        },
        std::vector
        {
        //      e,  A,  A,  C,  C,  G,  G,  T,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*A*/ N  ,DUL,DUL,N  ,N  ,N  ,N  ,N  ,N  ,N  ,DUL,DUL,N  ,N  ,N  ,N  ,N  ,N  ,
        /*C*/ N  ,N  ,N  ,DUL,DUL,N  ,N  ,N  ,N  ,N  ,N  ,N  ,DUl,DUL,N  ,N  ,N  ,N  ,
        /*G*/ N  ,N  ,N  ,N  ,DUL,DUL,DUL,N  ,N  ,N  ,N  ,N  ,N  ,DUl,DUL,DUL,N  ,N  ,
        /*T*/ N  ,N  ,N  ,N  ,N  ,N  ,DUL,DUL,DUL,DUl,N  ,N  ,N  ,N  ,N  ,DUl,DUL,DUL,
        /*C*/ N  ,N  ,N  ,DuL,DuL,N  ,N  ,N  ,DUl,N  ,N  ,N  ,Dul,DuL,N  ,N  ,N  ,DUl,
        /*T*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,DuL,DuL,DuL,N  ,N  ,N  ,N  ,N  ,N  ,Dul,DuL,
        /*A*/ N  ,DUL,DUL,N  ,N  ,N  ,N  ,N  ,N  ,N  ,DUL,DUL,N  ,N  ,N  ,N  ,N  ,N  ,
        /*C*/ N  ,N  ,N  ,DuL,DuL,N  ,N  ,N  ,N  ,N  ,Ul ,DUl,DuL,DuL,N  ,N  ,N  ,N  ,
        /*G*/ N  ,N  ,N  ,N  ,DUL,DuL,DUL,N  ,N  ,N  ,N  ,N  ,DUl,DUl,DuL,DUL,N  ,N  ,
        /*T*/ N  ,N  ,N  ,N  ,N  ,N  ,DUL,DuL,DuL,Dul,N  ,N  ,N  ,N  ,N  ,DUl,DuL,DuL,
        /*A*/ N  ,DuL,DuL,N  ,N  ,N  ,N  ,N  ,DUL,N  ,Dul,DuL,N  ,N  ,N  ,N  ,N  ,DUl
        }
    };
}();

// The same alignment with sequences swapped.
static auto dna4_02 = []()
{
    using seqan3::detail::column_index_type;
    using seqan3::detail::row_index_type;

    return alignment_fixture
    {
        "ACGTCTACGTA"_dna4,
        "AACCGGTTTAACCGGTT"_dna4,
        align_config | seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                    seqan3::mismatch_score{-5}}},
        11,
        "GTCTA",
        "GTTTA",
        seqan3::alignment_coordinate{column_index_type{2u}, row_index_type{5u}},
        seqan3::alignment_coordinate{column_index_type{7u}, row_index_type{10u}},
        std::vector
        {
        //     e, A, C, G, T, C, T, A, C, G, T, A
        /*e*/ 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
        /*A*/ 0 ,4 ,0 ,0 ,0 ,0 ,0 ,4 ,0 ,0 ,0 ,4 ,
        /*A*/ 0 ,4 ,0 ,0 ,0 ,0 ,0 ,4 ,0 ,0 ,0 ,4 ,
        /*C*/ 0 ,0 ,8 ,0 ,0 ,4 ,0 ,0 ,8 ,0 ,0 ,0 ,
        /*C*/ 0 ,0 ,4 ,3 ,0 ,4 ,0 ,0 ,4 ,3 ,0 ,0 ,
        /*G*/ 0 ,0 ,0 ,8 ,0 ,0 ,0 ,0 ,0 ,8 ,0 ,0 ,
        /*G*/ 0 ,0 ,0 ,4 ,3 ,0 ,0 ,0 ,0 ,4 ,3 ,0 ,
        /*T*/ 0 ,0 ,0 ,0 ,8 ,0 ,4 ,0 ,0 ,0 ,8 ,0 ,
        /*T*/ 0 ,0 ,0 ,0 ,4 ,3 ,4 ,0 ,0 ,0 ,4 ,3 ,
        /*T*/ 0 ,0 ,0 ,0 ,4 ,0 ,7 ,0 ,0 ,0 ,4 ,0 ,
        /*A*/ 0 ,4 ,0 ,0 ,0 ,0 ,0 ,11,0 ,0 ,0 ,8 ,
        /*A*/ 0 ,4 ,0 ,0 ,0 ,0 ,0 ,4 ,6 ,0 ,0 ,4 ,
        /*C*/ 0 ,0 ,8 ,0 ,0 ,4 ,0 ,0 ,8 ,1 ,0 ,0 ,
        /*C*/ 0 ,0 ,4 ,3 ,0 ,4 ,0 ,0 ,4 ,3 ,0 ,0 ,
        /*G*/ 0 ,0 ,0 ,8 ,0 ,0 ,0 ,0 ,0 ,8 ,0 ,0 ,
        /*G*/ 0 ,0 ,0 ,4 ,3 ,0 ,0 ,0 ,0 ,4 ,3 ,0 ,
        /*T*/ 0 ,0 ,0 ,0 ,8 ,0 ,4 ,0 ,0 ,0 ,8 ,0 ,
        /*T*/ 0 ,0 ,0 ,0 ,4 ,3 ,4 ,0 ,0 ,0 ,4 ,3
        },
        std::vector
        {
        //      e,  A,  C,  G,  T,  C,  T,  A,  C,  G,  T,  A
        /*e*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*A*/ N  ,DUL,N  ,N  ,N  ,N  ,N  ,DUL,N  ,N  ,N  ,DUl,
        /*A*/ N  ,DUL,N  ,N  ,N  ,N  ,N  ,DUL,N  ,N  ,N  ,DUl,
        /*C*/ N  ,N  ,DUL,N  ,N  ,DUl,N  ,N  ,DUl,N  ,N  ,N  ,
        /*C*/ N  ,N  ,DUL,DUL,N  ,DUl,N  ,N  ,DUl,DUL,N  ,N  ,
        /*G*/ N  ,N  ,N  ,DUL,N  ,N  ,N  ,N  ,N  ,DUl,N  ,N  ,
        /*G*/ N  ,N  ,N  ,DUL,DUL,N  ,N  ,N  ,N  ,DUL,DUL,N  ,
        /*T*/ N  ,N  ,N  ,N  ,DUL,N  ,DUl,N  ,N  ,N  ,DUl,N  ,
        /*T*/ N  ,N  ,N  ,N  ,DUL,DuL,DUl,N  ,N  ,N  ,DUl,DUL,
        /*T*/ N  ,N  ,N  ,N  ,DuL,N  ,DUl,N  ,N  ,N  ,Dul,N  ,
        /*A*/ N  ,DUL,N  ,N  ,N  ,N  ,N  ,DUL,L  ,N  ,N  ,Dul,
        /*A*/ N  ,DUL,N  ,N  ,N  ,N  ,N  ,DUL,DuL,N  ,N  ,DUl,
        /*C*/ N  ,N  ,DuL,N  ,N  ,Dul,N  ,N  ,DUl,DuL,N  ,N  ,
        /*C*/ N  ,N  ,DUL,DuL,N  ,DUl,N  ,N  ,DUl,DuL,N  ,N  ,
        /*G*/ N  ,N  ,N  ,DUL,N  ,N  ,N  ,N  ,N  ,DUl,N  ,N  ,
        /*G*/ N  ,N  ,N  ,DUL,DuL,N  ,N  ,N  ,N  ,DUL,DuL,N  ,
        /*T*/ N  ,N  ,N  ,N  ,DUL,N  ,Dul,N  ,N  ,N  ,DUl,N  ,
        /*T*/ N  ,N  ,N  ,N  ,DUL,DuL,DUl,N  ,N  ,N  ,DUl,DuL
        }
    };
}();

// Local alignment starting in the first row. Verifies that free end gaps are performed correctly.
static auto dna4_03 = []()
{
    using seqan3::detail::column_index_type;
    using seqan3::detail::row_index_type;

    return alignment_fixture
    {
        "ataagcgtctcg"_dna4,
        "tcatagagttgc"_dna4,
        seqan3::align_cfg::mode{seqan3::local_alignment}
            | seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1}, seqan3::gap_open_score{-1}}}
            | seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{seqan3::match_score{2},
                                                                           seqan3::mismatch_score{-1}}},
        9,
        "ATAAGCGT",
        "AT-AGAGT",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{2u}},
        seqan3::alignment_coordinate{column_index_type{8u}, row_index_type{9u}},
        std::vector
        {
        //    e,A,T,A,A,G,C,G,T,C,T,C,G
        /*e*/ 0,0,0,0,0,0,0,0,0,0,0,0,0,
        /*T*/ 0,0,2,0,0,0,0,0,2,0,2,0,0,
        /*C*/ 0,0,0,1,0,0,2,0,0,4,2,4,2,
        /*A*/ 0,2,0,2,3,1,0,1,0,2,3,2,3,
        /*T*/ 0,0,4,2,1,2,0,0,3,1,4,2,1,
        /*A*/ 0,2,2,6,4,3,2,1,1,2,2,3,1,
        /*G*/ 0,0,1,4,5,6,4,4,2,1,1,1,5,
        /*A*/ 0,2,0,3,6,4,5,3,3,1,0,0,3,
        /*G*/ 0,0,1,2,4,8,6,7,5,4,3,2,2,
        /*T*/ 0,0,2,1,3,6,7,5,9,7,6,5,4,
        /*T*/ 0,0,2,1,2,5,5,6,7,8,9,7,6,
        /*G*/ 0,0,0,1,1,4,4,7,6,6,7,8,9,
        /*C*/ 0,0,0,0,0,3,6,5,6,8,6,9,7
        },
        std::vector
        {
        //      e,  A,  T,  A,  A,  G,  C,  G,  T,  C,  T,  C,  G
        /*e*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*T*/ N  ,N  ,DUL,L  ,N  ,N  ,N  ,N  ,DUL,L  ,DUl,L  ,N  ,
        /*C*/ N  ,N  ,UL ,DUL,N  ,N  ,DUL,L  ,Ul ,DUl,L  ,DUl,L  ,
        /*A*/ N  ,DUL,L  ,DUl,DUL,L  ,Ul ,DUl,N  ,Ul ,DUL,UL ,DUl,
        /*T*/ N  ,UL ,DuL,L  ,DUl,DUl,DuL,N  ,Dul,uL ,DUl,DuL,DUl,
        /*A*/ N  ,DuL,UL ,DUL,DuL,l  ,l  ,l  ,Ul ,Dul,UL ,DuL,DuL,
        /*G*/ N  ,UL ,DuL,UL ,DUL,DUL,L  ,DUl,l  ,l  ,Dul,DUl,DuL,
        /*A*/ N  ,DuL,uL ,Dul,DUL,DUL,DUl,DUL,DUl,Dul,Dul,Dul,Ul ,
        /*G*/ N  ,UL ,DuL,uL ,UL ,DuL,L  ,Dul,L  ,l  ,l  ,l  ,Dul,
        /*T*/ N  ,N  ,DUL,uL ,ul ,UL ,DUL,DUL,DUl,L  ,DUl,l  ,l  ,
        /*T*/ N  ,N  ,DUL,DuL,ul ,uL ,DUL,DuL,DUL,DUL,DUL,L  ,l  ,
        /*G*/ N  ,N  ,UL ,DuL,uL ,DuL,DuL,DUL,uL ,DUl,DUL,DUL,DUL,
        /*C*/ N  ,N  ,N  ,N  ,DuL,uL ,DuL,UL ,Dul,DuL,uL ,DUl,DUL
        }
    };
}();

// Only mismatches, so an empty alignment is found (score 0).
static auto dna4_04 = []()
{
    using seqan3::detail::column_index_type;
    using seqan3::detail::row_index_type;

    return alignment_fixture
    {
        "AAAAAA"_dna4,
        "CCCCCC"_dna4,
        align_config | seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                    seqan3::mismatch_score{-5}}},
        0,
        "",
        "",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        std::vector
        {
        //    e,A,A,A,A,A,A
        /*e*/ 0,0,0,0,0,0,0,
        /*C*/ 0,0,0,0,0,0,0,
        /*C*/ 0,0,0,0,0,0,0,
        /*C*/ 0,0,0,0,0,0,0,
        /*C*/ 0,0,0,0,0,0,0,
        /*C*/ 0,0,0,0,0,0,0,
        /*C*/ 0,0,0,0,0,0,0
        },
        std::vector
        {
        //    e,A,A,A,A,A,A
        /*e*/ N,N,N,N,N,N,N,
        /*C*/ N,N,N,N,N,N,N,
        /*C*/ N,N,N,N,N,N,N,
        /*C*/ N,N,N,N,N,N,N,
        /*C*/ N,N,N,N,N,N,N,
        /*C*/ N,N,N,N,N,N,N,
        /*C*/ N,N,N,N,N,N,N
        }
    };
}();

// Local alignment in the begin and end of sequences.
static auto dna4_05 = []()
{
    using seqan3::detail::column_index_type;
    using seqan3::detail::row_index_type;

    return alignment_fixture
    {
        "AAAAAATCCCCCC"_dna4,
        "CCCCCCTAAAAAA"_dna4,
        align_config | seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                    seqan3::mismatch_score{-5}}},
        24,
        "AAAAAA",
        "AAAAAA",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{7u}},
        seqan3::alignment_coordinate{column_index_type{6u}, row_index_type{13u}},
        std::vector
        {
        //     e, A, A, A, A, A, A, T, C, C, C, C, C, C
        /*e*/ 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
        /*C*/ 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4 ,4 ,4 ,4 ,4 ,4 ,
        /*C*/ 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4 ,8 ,8 ,8 ,8 ,8 ,
        /*C*/ 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4 ,8 ,12,12,12,12,
        /*C*/ 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4 ,8 ,12,16,16,16,
        /*C*/ 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4 ,8 ,12,16,20,20,
        /*C*/ 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4 ,8 ,12,16,20,24,
        /*T*/ 0 ,0 ,0 ,0 ,0 ,0 ,0 ,4 ,0 ,0 ,3 ,7 ,11,15,
        /*A*/ 0 ,4 ,4 ,4 ,4 ,4 ,4 ,0 ,0 ,0 ,0 ,4 ,8 ,12,
        /*A*/ 0 ,4 ,8 ,8 ,8 ,8 ,8 ,0 ,0 ,0 ,0 ,3 ,7 ,11,
        /*A*/ 0 ,4 ,8 ,12,12,12,12,3 ,0 ,0 ,0 ,2 ,6 ,10,
        /*A*/ 0 ,4 ,8 ,12,16,16,16,7 ,4 ,3 ,2 ,1 ,5 ,9 ,
        /*A*/ 0 ,4 ,8 ,12,16,20,20,11,8 ,7 ,6 ,5 ,4 ,8 ,
        /*A*/ 0 ,4 ,8 ,12,16,20,24,15,12,11,10,9 ,8 ,7
        },
        std::vector
        {
        //      e,  A,  A,  A,  A,  A,  A,  T,  C,  C,  C,  C,  C,  C
        /*e*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*C*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,DUL,DUL,DUL,DUL,DUL,DUL,
        /*C*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,DUL,DUL,DUL,DUL,DUL,DUL,
        /*C*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,DUL,DUL,DUL,DUL,DUL,DUL,
        /*C*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,DUL,DUL,DUL,DUL,DUL,DUL,
        /*C*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,DUL,DUL,DUL,DUL,DUL,DUL,
        /*C*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,DUL,DUL,DUL,DUL,DUL,DUL,
        /*T*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,DUL,N  ,N  ,DUl,DUL,DUL,DUL,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,N  ,N  ,N  ,ul ,ul ,uL ,uL ,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,N  ,N  ,N  ,N  ,ul ,ul ,uL ,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DuL,l  ,N  ,N  ,ul ,ul ,ul ,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUL,l  ,l  ,l  ,ul ,ul ,ul ,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUL,l  ,l  ,l  ,l  ,ul ,ul ,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUL,l  ,l  ,l  ,l  ,l  ,ul
        }
    };
}();

// Local RNA alignment with a longer sequence of gaps.
static auto rna5_01 = []()
{
    using seqan3::detail::column_index_type;
    using seqan3::detail::row_index_type;

    return alignment_fixture
    {
        "AAAAAAUUUUNNUUUUCCCCCC"_rna5,
        "AAAAAACCCCCC"_rna5,
        align_config | seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                    seqan3::mismatch_score{-5}}},
        28,
        "AAAAAAUUUUNNUUUUCCCCCC",
        "AAAAAA----------CCCCCC",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        seqan3::alignment_coordinate{column_index_type{22u}, row_index_type{12u}},
        std::vector
        {
        //     e, A, A, A, A, A, A, U, U, U, U, N, N, U, U, U, U, C, C, C, C, C, C
        /*e*/ 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
        /*A*/ 0 ,4 ,4 ,4 ,4 ,4 ,4 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
        /*A*/ 0 ,4 ,8 ,8 ,8 ,8 ,8 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
        /*A*/ 0 ,4 ,8 ,12,12,12,12,3 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
        /*A*/ 0 ,4 ,8 ,12,16,16,16,7 ,4 ,3 ,2 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
        /*A*/ 0 ,4 ,8 ,12,16,20,20,11,8 ,7 ,6 ,5 ,4 ,3 ,2 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
        /*A*/ 0 ,4 ,8 ,12,16,20,24,15,12,11,10,9 ,8 ,7 ,6 ,5 ,4 ,3 ,2 ,1 ,0 ,0 ,0 ,
        /*C*/ 0 ,0 ,0 ,3 ,7 ,11,15,19,10,7 ,6 ,5 ,4 ,3 ,2 ,1 ,0 ,8 ,7 ,6 ,5 ,4 ,4 ,
        /*C*/ 0 ,0 ,0 ,0 ,4 ,8 ,12,10,14,5 ,2 ,1 ,0 ,0 ,0 ,0 ,0 ,4 ,12,11,10,9 ,8 ,
        /*C*/ 0 ,0 ,0 ,0 ,3 ,7 ,11,7 ,5 ,9 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4 ,8 ,16,15,14,13,
        /*C*/ 0 ,0 ,0 ,0 ,2 ,6 ,10,6 ,2 ,0 ,4 ,0 ,0 ,0 ,0 ,0 ,0 ,4 ,8 ,12,20,19,18,
        /*C*/ 0 ,0 ,0 ,0 ,1 ,5 ,9 ,5 ,1 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4 ,8 ,12,16,24,23,
        /*C*/ 0 ,0 ,0 ,0 ,0 ,4 ,8 ,4 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,4 ,8 ,12,16,20,28
        },
        std::vector
        {
        //      e,  A,  A,  A,  A,  A,  A,  U,  U,  U,  U,  N,  N,  U,  U,  U,  U,  C,  C,  C,  C,  C,  C
        /*e*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUL,l  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUL,l  ,l  ,l  ,l  ,l  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUL,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUL,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,N  ,N  ,
        /*C*/ N  ,N  ,N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,
        /*C*/ N  ,N  ,N  ,uL ,uL ,uL ,uL ,DUL,Dul,DuL,Dul,Dul,Dul,N  ,N  ,N  ,N  ,DUl,DUl,DUL,DUl,DUl,DUl,
        /*C*/ N  ,N  ,N  ,N  ,uL ,uL ,uL ,DuL,DUl,Dul,DuL,N  ,N  ,N  ,N  ,N  ,N  ,Dul,DUL,DUL,DUL,DUl,DUl,
        /*C*/ N  ,N  ,N  ,N  ,uL ,uL ,uL ,DuL,Dul,DUl,Dul,N  ,N  ,N  ,N  ,N  ,N  ,Dul,DuL,DUL,DUL,DUL,DUl,
        /*C*/ N  ,N  ,N  ,N  ,uL ,uL ,uL ,DuL,Dul,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,DuL,DuL,DuL,DUL,DUL,DUL,
        /*C*/ N  ,N  ,N  ,N  ,uL ,uL ,uL ,DuL,Dul,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,DuL,DuL,DuL,DuL,DUL,DUL
        }
    };
}();

// Local alignment for proteins (amino acid sequence) with BLOSUM62 score.
static auto aa27_01 = []()
{
    using seqan3::detail::column_index_type;
    using seqan3::detail::row_index_type;

    return alignment_fixture
    {
        "ALIGATOR"_aa27,
        "GALORA"_aa27,
        align_config | seqan3::align_cfg::scoring{aminoacid_scoring_scheme{aminoacid_similarity_matrix::BLOSUM62}},
        13,
        "GATOR",
        "GALOR",
        seqan3::alignment_coordinate{column_index_type{3u}, row_index_type{0u}},
        seqan3::alignment_coordinate{column_index_type{8u}, row_index_type{5u}},
        std::vector
        {
        //     e, A, L, I, G, A, T, O, R
        /*e*/ 0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,0 ,
        /*G*/ 0 ,0 ,0 ,0 ,6 ,0 ,0 ,0 ,0 ,
        /*A*/ 0 ,4 ,0 ,0 ,0 ,10,0 ,0 ,0 ,
        /*L*/ 0 ,0 ,8 ,2 ,0 ,0 ,9 ,0 ,0 ,
        /*O*/ 0 ,0 ,0 ,7 ,1 ,0 ,0 ,8 ,0 ,
        /*R*/ 0 ,0 ,0 ,0 ,5 ,0 ,0 ,0 ,13,
        /*A*/ 0 ,4 ,0 ,0 ,0 ,9 ,0 ,0 ,2
        },
        std::vector
        {
        //      e,  A,  L,  I,  G,  A,  T,  O,  R
        /*e*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*G*/ N  ,DUL,N  ,N  ,DUL,DUL,N  ,N  ,N  ,
        /*A*/ N  ,DUL,N  ,N  ,DUl,DUl,DUL,DUl,N  ,
        /*L*/ N  ,N  ,DUL,DUL,N  ,N  ,DUl,N  ,N  ,
        /*O*/ N  ,DuL,N  ,DUL,DuL,Dul,DUl,DUl,N  ,
        /*R*/ N  ,N  ,N  ,N  ,DuL,DuL,N  ,N  ,DUl,
        /*A*/ N  ,DuL,N  ,N  ,DUl,Dul,DuL,Dul,Ul
        }
    };
}();

// Local alignment with empty sequence.
static auto aa27_02 = []()
{
    using seqan3::detail::column_index_type;
    using seqan3::detail::row_index_type;

    return alignment_fixture
    {
        "ALIGATOR"_aa27,
        ""_aa27,
        align_config | seqan3::align_cfg::scoring{aminoacid_scoring_scheme{aminoacid_similarity_matrix::BLOSUM62}},
        0,
        "",
        "",
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        seqan3::alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        std::vector
        {
        //e,A,L,I,G,A,T,O,R
          0,0,0,0,0,0,0,0,0
        },
        std::vector
        {
        //e,A,L,I,G,A,T,O,R
          N,N,N,N,N,N,N,N,N
        }
    };
}();

} // namespace seqan3::test::alignment::fixture::local::affine::unbanded
