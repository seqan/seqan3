// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <vector>

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_mode.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>

#include "alignment_fixture.hpp"

namespace seqan3::test::alignment::fixture::local::affine::banded
{

inline constexpr auto align_config = align_cfg::mode{local_alignment} |
                                     align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}};

// Local alignment with mismatch.
static auto dna4_01 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        // score: 11 (4 matches, 1 mismatch)
        // alignment:
        // GTTTA
        // || ||
        // GTCTA
        "AACCGGTTTAACCGGTT"_dna4,
        "ACGTCTACGTA"_dna4,
        align_config
            | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}}
            | align_cfg::band{static_band{lower_bound{-2}, upper_bound{5}}},
        11,
        "GTTTA",
        "GTCTA",
        alignment_coordinate{column_index_type{5u}, row_index_type{2u}},
        alignment_coordinate{column_index_type{10u}, row_index_type{7u}},
        std::vector<std::optional<int32_t>>
        {
        //      e,  A,  A,  C,  C,  G,  G,  T,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/ 0  ,0  ,0  ,0  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ 0  ,4  ,4  ,0  ,0  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ 0  ,0  ,0  ,8  ,4  ,0  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*G*/ INF,0  ,0  ,0  ,3  ,8  ,4  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ INF,INF,0  ,0  ,0  ,0  ,3  ,8  ,4  ,4  ,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ INF,INF,INF,4  ,4  ,0  ,0  ,0  ,3  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ INF,INF,INF,INF,0  ,0  ,0  ,4  ,4  ,7  ,0  ,0  ,INF,INF,INF,INF,INF,INF,
        /*A*/ INF,INF,INF,INF,INF,0  ,0  ,0  ,0  ,0  ,11 ,4  ,0  ,INF,INF,INF,INF,INF,
        /*C*/ INF,INF,INF,INF,INF,INF,0  ,0  ,0  ,0  ,0  ,6  ,8  ,4  ,INF,INF,INF,INF,
        /*G*/ INF,INF,INF,INF,INF,INF,INF,0  ,0  ,0  ,0  ,0  ,1  ,3  ,8  ,INF,INF,INF,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,INF,4  ,4  ,0  ,0  ,0  ,0  ,0  ,3  ,INF,INF,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,0  ,8  ,4  ,0  ,0  ,0  ,0  ,0  ,INF
        },
        std::vector<std::optional<detail::trace_directions>>
        {
        //      e,  A,  A,  C,  C,  G,  G,  T,  T,  T,  A,  A,  C,  C,  G,  G,  T,  T
        /*e*/ N  ,N  ,N  ,N  ,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ N  ,DUL,DUL,N  ,N  ,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ N  ,N  ,N  ,DUL,DUL,N  ,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*G*/ INF,N  ,N  ,N  ,DUL,DUL,DUL,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ INF,INF,N  ,N  ,N  ,N  ,DUL,DUL,DUL,D  ,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ INF,INF,INF,Du ,DuL,N  ,N  ,N  ,DUl,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ INF,INF,INF,INF,N  ,N  ,N  ,DuL,DuL,DuL,N  ,N  ,INF,INF,INF,INF,INF,INF,
        /*A*/ INF,INF,INF,INF,INF,N  ,N  ,N  ,N  ,N  ,DUL,DUL,N  ,INF,INF,INF,INF,INF,
        /*C*/ INF,INF,INF,INF,INF,INF,N  ,N  ,N  ,N  ,UL ,DUL,DUL,D  ,INF,INF,INF,INF,
        /*G*/ INF,INF,INF,INF,INF,INF,INF,N  ,N  ,N  ,N  ,N  ,DUL,DUL,D  ,INF,INF,INF,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,INF,Du ,DuL,N  ,N  ,N  ,N  ,N  ,D  ,INF,INF,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,N  ,DuL,DuL,N  ,N  ,N  ,N  ,N  ,INF
        }
    };
}();

// The same alignment with sequences swapped. The asymmetric band leads to a worse result than above.
static auto dna4_02 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        "ACGTCTACGTA"_dna4,
        "AACCGGTTTAACCGGTT"_dna4,
        align_config
            | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}}
            | align_cfg::band{static_band{lower_bound{-2}, upper_bound{5}}},
        8,
        "AC",
        "AC",
        alignment_coordinate{column_index_type{0u}, row_index_type{1u}},
        alignment_coordinate{column_index_type{2u}, row_index_type{3u}},
        std::vector<std::optional<int32_t>>
        {
        //      e,  A,  C,  G,  T,  C,  T,  A,  C,  G,  T,  A
        /*e*/ 0  ,0  ,0  ,0  ,0  ,0  ,INF,INF,INF,INF,INF,INF,
        /*A*/ 0  ,4  ,0  ,0  ,0  ,0  ,0  ,INF,INF,INF,INF,INF,
        /*A*/ 0  ,4  ,0  ,0  ,0  ,0  ,0  ,4  ,INF,INF,INF,INF,
        /*C*/ INF,0  ,8  ,0  ,0  ,4  ,0  ,0  ,8  ,INF,INF,INF,
        /*C*/ INF,INF,4  ,3  ,0  ,4  ,0  ,0  ,4  ,3  ,INF,INF,
        /*G*/ INF,INF,INF,8  ,0  ,0  ,0  ,0  ,0  ,8  ,0  ,INF,
        /*G*/ INF,INF,INF,INF,3  ,0  ,0  ,0  ,0  ,4  ,3  ,0  ,
        /*T*/ INF,INF,INF,INF,INF,0  ,4  ,0  ,0  ,0  ,8  ,0  ,
        /*T*/ INF,INF,INF,INF,INF,INF,4  ,0  ,0  ,0  ,4  ,3  ,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,0  ,0  ,0  ,4  ,0  ,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,0  ,0  ,0  ,8  ,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,0  ,0  ,4  ,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,0  ,0  ,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,0  ,
        /*G*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*G*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF
        },
        std::vector<std::optional<detail::trace_directions>>
        {
        //      e,  A,  C,  G,  T,  C,  T,  A,  C,  G,  T,  A
        /*e*/ N  ,N  ,N  ,N  ,N  ,N  ,INF,INF,INF,INF,INF,INF,
        /*A*/ N  ,DUL,N  ,N  ,N  ,N  ,N  ,INF,INF,INF,INF,INF,
        /*A*/ N  ,DUL,N  ,N  ,N  ,N  ,N  ,D  ,INF,INF,INF,INF,
        /*C*/ INF,N  ,DUL,N  ,N  ,DUl,N  ,N  ,D  ,INF,INF,INF,
        /*C*/ INF,INF,DU ,DUL,N  ,DUl,N  ,N  ,DUl,D  ,INF,INF,
        /*G*/ INF,INF,INF,DU ,N  ,N  ,N  ,N  ,N  ,DUl,N  ,INF,
        /*G*/ INF,INF,INF,INF,DU ,N  ,N  ,N  ,N  ,DUL,DUL,N  ,
        /*T*/ INF,INF,INF,INF,INF,N  ,DUL,N  ,N  ,N  ,DUl,N  ,
        /*T*/ INF,INF,INF,INF,INF,INF,DU ,N  ,N  ,N  ,DUl,DUL,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,N  ,N  ,N  ,DuL,N  ,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,N  ,N  ,N  ,DuL,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,N  ,N  ,DUL,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,N  ,N  ,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,N  ,
        /*G*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*G*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF
        }
    };
}();

// Local alignment with zero bandwidth. Does not allow any gaps.
static auto dna4_03 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        "ataagcgtctcg"_dna4,
        "ctcagagttgc"_dna4,
        align_cfg::mode{local_alignment}
            | align_cfg::gap{gap_scheme{gap_score{0}, gap_open_score{0}}}
            | align_cfg::scoring{nucleotide_scoring_scheme{match_score{2}, mismatch_score{-1}}}
            | align_cfg::band{static_band{lower_bound{0}, upper_bound{0}}},
        8,
        "TAAGCGT",
        "TCAGAGT",
        alignment_coordinate{column_index_type{1u}, row_index_type{1u}},
        alignment_coordinate{column_index_type{8u}, row_index_type{8u}},
        std::vector<std::optional<int32_t>>
        {
        //      e,  A,  T,  A,  A,  G,  C,  G,  T,  C,  T,  C,  G
        /*e*/ 0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ INF,0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ INF,INF,2  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ INF,INF,INF,1  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ INF,INF,INF,INF,3  ,INF,INF,INF,INF,INF,INF,INF,INF,
        /*G*/ INF,INF,INF,INF,INF,5  ,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ INF,INF,INF,INF,INF,INF,4  ,INF,INF,INF,INF,INF,INF,
        /*G*/ INF,INF,INF,INF,INF,INF,INF,6  ,INF,INF,INF,INF,INF,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,INF,8  ,INF,INF,INF,INF,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,7  ,INF,INF,INF,
        /*G*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,6  ,INF,INF,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,8  ,INF
        },
        std::vector<std::optional<detail::trace_directions>>
        {
        //      e,  A,  T,  A,  A,  G,  C,  G,  T,  C,  T,  C,  G
        /*e*/ N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ INF,N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ INF,INF,D  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ INF,INF,INF,D  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ INF,INF,INF,INF,D  ,INF,INF,INF,INF,INF,INF,INF,INF,
        /*G*/ INF,INF,INF,INF,INF,D  ,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ INF,INF,INF,INF,INF,INF,D  ,INF,INF,INF,INF,INF,INF,
        /*G*/ INF,INF,INF,INF,INF,INF,INF,D  ,INF,INF,INF,INF,INF,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,INF,D  ,INF,INF,INF,INF,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,D  ,INF,INF,INF,
        /*G*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,D  ,INF,INF,
        /*C*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,D  ,INF
        }
    };
}();

// Only mismatches, so an empty alignment is found (score 0).
static auto dna4_04 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        "AAAAAA"_dna4,
        "CCCCCC"_dna4,
        align_config
            | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}}
            | align_cfg::band{static_band{lower_bound{-2}, upper_bound{2}}},
        0,
        "",
        "",
        alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        std::vector<std::optional<int32_t>>
        {
        //      e,  A,  A,  A,  A,  A,  A
        /*e*/ 0  ,0  ,0  ,INF,INF,INF,INF,
        /*C*/ 0  ,0  ,0  ,0  ,INF,INF,INF,
        /*C*/ 0  ,0  ,0  ,0  ,0  ,INF,INF,
        /*C*/ INF,0  ,0  ,0  ,0  ,0  ,INF,
        /*C*/ INF,INF,0  ,0  ,0  ,0  ,0  ,
        /*C*/ INF,INF,INF,0  ,0  ,0  ,0  ,
        /*C*/ INF,INF,INF,INF,0  ,0  ,0
        },
        std::vector<std::optional<detail::trace_directions>>
        {
        //      e,  A,  A,  A,  A,  A,  A
        /*e*/ N  ,N  ,N  ,INF,INF,INF,INF,
        /*C*/ N  ,N  ,N  ,N  ,INF,INF,INF,
        /*C*/ N  ,N  ,N  ,N  ,N  ,INF,INF,
        /*C*/ INF,N  ,N  ,N  ,N  ,N  ,INF,
        /*C*/ INF,INF,N  ,N  ,N  ,N  ,N  ,
        /*C*/ INF,INF,INF,N  ,N  ,N  ,N  ,
        /*C*/ INF,INF,INF,INF,N  ,N  ,N
        }
    };
}();

// Local alignment in the begin and end of sequences. The band covers the lower diagonal matrix.
static auto dna4_05 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        "AAAAAATCCCCCC"_dna4,
        "CCCCCCTAAAAAA"_dna4,
        align_config
            | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}}
            | align_cfg::band{static_band{lower_bound{-100}, upper_bound{0}}},
        24,
        "AAAAAA",
        "AAAAAA",
        alignment_coordinate{column_index_type{0u}, row_index_type{7u}},
        alignment_coordinate{column_index_type{6u}, row_index_type{13u}},
        std::vector<std::optional<int32_t>>
        {
        //      e,  A,  A,  A,  A,  A,  A,  T,  C,  C,  C,  C,  C,  C
        /*e*/ 0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ 0  ,0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ 0  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ 0  ,0  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ 0  ,0  ,0  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ 0  ,0  ,0  ,0  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ 0  ,0  ,0  ,0  ,0  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ 0  ,0  ,0  ,0  ,0  ,0  ,0  ,4  ,INF,INF,INF,INF,INF,INF,
        /*A*/ 0  ,4  ,4  ,4  ,4  ,4  ,4  ,0  ,0  ,INF,INF,INF,INF,INF,
        /*A*/ 0  ,4  ,8  ,8  ,8  ,8  ,8  ,0  ,0  ,0  ,INF,INF,INF,INF,
        /*A*/ 0  ,4  ,8  ,12 ,12 ,12 ,12 ,3  ,0  ,0  ,0  ,INF,INF,INF,
        /*A*/ 0  ,4  ,8  ,12 ,16 ,16 ,16 ,7  ,4  ,3  ,2  ,1  ,INF,INF,
        /*A*/ 0  ,4  ,8  ,12 ,16 ,20 ,20 ,11 ,8  ,7  ,6  ,5  ,4  ,INF,
        /*A*/ 0  ,4  ,8  ,12 ,16 ,20 ,24 ,15 ,12 ,11 ,10 ,9  ,8  ,7
        },
        std::vector<std::optional<detail::trace_directions>>
        {
        //      e,  A,  A,  A,  A,  A,  A,  T,  C,  C,  C,  C,  C,  C
        /*e*/ N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ N  ,N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ N  ,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ N  ,N  ,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ N  ,N  ,N  ,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ N  ,N  ,N  ,N  ,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,INF,
        /*C*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,
        /*T*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,D  ,INF,INF,INF,INF,INF,INF,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,N  ,N  ,INF,INF,INF,INF,INF,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,N  ,N  ,N  ,INF,INF,INF,INF,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DuL,l  ,N  ,N  ,INF,INF,INF,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUL,l  ,l  ,l  ,l  ,INF,INF,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUL,l  ,l  ,l  ,l  ,l  ,INF,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUL,l  ,l  ,l  ,l  ,l  ,l
        }
    };
}();

// Local alignment in the begin and end of sequences. The band cover the upper diagonal matrix and
// enforces to align the C's instead of the A's.
static auto dna4_06 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        "AAAAAATCCCCCC"_dna4,
        "CCCCCCTAAAAAA"_dna4,
        align_config
            | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}}
            | align_cfg::band{static_band{lower_bound{0}, upper_bound{100}}},
        24,
        "CCCCCC",
        "CCCCCC",
        alignment_coordinate{column_index_type{7u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{13u}, row_index_type{6u}},
        std::vector<std::optional<int32_t>>
        {
        //      e,  A,  A,  A,  A,  A,  A,  T,  C,  C,  C,  C,  C,  C
        /*e*/ 0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,
        /*C*/ INF,0  ,0  ,0  ,0  ,0  ,0  ,0  ,4  ,4  ,4  ,4  ,4  ,4  ,
        /*C*/ INF,INF,0  ,0  ,0  ,0  ,0  ,0  ,4  ,8  ,8  ,8  ,8  ,8  ,
        /*C*/ INF,INF,INF,0  ,0  ,0  ,0  ,0  ,4  ,8  ,12 ,12 ,12 ,12 ,
        /*C*/ INF,INF,INF,INF,0  ,0  ,0  ,0  ,4  ,8  ,12 ,16 ,16 ,16 ,
        /*C*/ INF,INF,INF,INF,INF,0  ,0  ,0  ,4  ,8  ,12 ,16 ,20 ,20 ,
        /*C*/ INF,INF,INF,INF,INF,INF,0  ,0  ,4  ,8  ,12 ,16 ,20 ,24 ,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,4  ,0  ,0  ,3  ,7  ,11 ,15 ,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,0  ,0  ,0  ,4  ,8  ,12 ,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,0  ,0  ,3  ,7  ,11 ,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,0  ,2  ,6  ,10 ,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,1  ,5  ,9  ,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,4  ,8  ,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,7
        },
        std::vector<std::optional<detail::trace_directions>>
        {
        //      e,  A,  A,  A,  A,  A,  A,  T,  C,  C,  C,  C,  C,  C
        /*e*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,
        /*C*/ INF,N  ,N  ,N  ,N  ,N  ,N  ,N  ,DUL,DUL,DUL,DUL,DUL,DUL,
        /*C*/ INF,INF,N  ,N  ,N  ,N  ,N  ,N  ,DUL,DUL,DUL,DUL,DUL,DUL,
        /*C*/ INF,INF,INF,N  ,N  ,N  ,N  ,N  ,DUL,DUL,DUL,DUL,DUL,DUL,
        /*C*/ INF,INF,INF,INF,N  ,N  ,N  ,N  ,DUL,DUL,DUL,DUL,DUL,DUL,
        /*C*/ INF,INF,INF,INF,INF,N  ,N  ,N  ,DUL,DUL,DUL,DUL,DUL,DUL,
        /*C*/ INF,INF,INF,INF,INF,INF,N  ,N  ,DUL,DUL,DUL,DUL,DUL,DUL,
        /*T*/ INF,INF,INF,INF,INF,INF,INF,DU ,N  ,N  ,DUl,DUL,DUL,DUL,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,N  ,N  ,uL ,uL ,uL ,uL ,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,N  ,N  ,uL ,uL ,uL ,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,N  ,uL ,uL ,uL ,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,u  ,uL ,uL ,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,u  ,uL ,
        /*A*/ INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,u
        }
    };
}();

// Local RNA alignment with a longer sequence of gaps. The alignment trace is located along the band boundary.
static auto rna5_01 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        "AAAAAAUUUUNNUUUUCCCCCC"_rna5,
        "AAAAAACCCCCC"_rna5,
        align_config
            | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}}
            | align_cfg::band{static_band{lower_bound{-10}, upper_bound{10}}},
        28,
        "AAAAAAUUUUNNUUUUCCCCCC",
        "AAAAAA----------CCCCCC",
        alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{22u}, row_index_type{12u}},
        std::vector<std::optional<int32_t>>
        {
        //      e,  A,  A,  A,  A,  A,  A,  U,  U,  U,  U,  N,  N,  U,  U,  U,  U,  C,  C,  C,  C,  C,  C
        /*e*/ 0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ 0  ,4  ,4  ,4  ,4  ,4  ,4  ,0  ,0  ,0  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ 0  ,4  ,8  ,8  ,8  ,8  ,8  ,0  ,0  ,0  ,0  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ 0  ,4  ,8  ,12 ,12 ,12 ,12 ,3  ,0  ,0  ,0  ,0  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ 0  ,4  ,8  ,12 ,16 ,16 ,16 ,7  ,4  ,3  ,2  ,1  ,0  ,0  ,0  ,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ 0  ,4  ,8  ,12 ,16 ,20 ,20 ,11 ,8  ,7  ,6  ,5  ,4  ,3  ,2  ,1  ,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ 0  ,4  ,8  ,12 ,16 ,20 ,24 ,15 ,12 ,11 ,10 ,9  ,8  ,7  ,6  ,5  ,4  ,INF,INF,INF,INF,INF,INF,
        /*C*/ 0  ,0  ,0  ,3  ,7  ,11 ,15 ,19 ,10 ,7  ,6  ,5  ,4  ,3  ,2  ,1  ,0  ,8  ,INF,INF,INF,INF,INF,
        /*C*/ 0  ,0  ,0  ,0  ,4  ,8  ,12 ,10 ,14 ,5  ,2  ,1  ,0  ,0  ,0  ,0  ,0  ,4  ,12 ,INF,INF,INF,INF,
        /*C*/ 0  ,0  ,0  ,0  ,3  ,7  ,11 ,7  ,5  ,9  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,4  ,8  ,16 ,INF,INF,INF,
        /*C*/ 0  ,0  ,0  ,0  ,2  ,6  ,10 ,6  ,2  ,0  ,4  ,0  ,0  ,0  ,0  ,0  ,0  ,4  ,8  ,12 ,20 ,INF,INF,
        /*C*/ INF,0  ,0  ,0  ,1  ,5  ,9  ,5  ,1  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,4  ,8  ,12 ,16 ,24 ,INF,
        /*C*/ INF,INF,0  ,0  ,0  ,4  ,8  ,4  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,0  ,4  ,8  ,12 ,16 ,20 ,28
        },
        std::vector<std::optional<detail::trace_directions>>
        {
        //      e,  A,  A,  A,  A,  A,  A,  U,  U,  U,  U,  N,  N,  U,  U,  U,  U,  C,  C,  C,  C,  C,  C
        /*e*/ N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,N  ,N  ,N  ,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,N  ,N  ,N  ,N  ,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUL,l  ,N  ,N  ,N  ,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUL,l  ,l  ,l  ,l  ,l  ,N  ,N  ,INF,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUL,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,INF,INF,INF,INF,INF,INF,INF,
        /*A*/ N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUL,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,INF,INF,INF,INF,INF,INF,
        /*C*/ N  ,N  ,N  ,DUL,DUL,DUL,DUL,DUL,DUL,DUl,DUl,DUl,DUl,DUl,DUl,DUl,DUl,D  ,INF,INF,INF,INF,INF,
        /*C*/ N  ,N  ,N  ,uL ,uL ,uL ,uL ,DUL,Dul,DuL,Dul,Dul,Dul,N  ,N  ,N  ,N  ,DUl,D  ,INF,INF,INF,INF,
        /*C*/ N  ,N  ,N  ,N  ,uL ,uL ,uL ,DuL,DUl,Dul,DuL,N  ,N  ,N  ,N  ,N  ,N  ,Dul,DUL,D  ,INF,INF,INF,
        /*C*/ N  ,N  ,N  ,N  ,uL ,uL ,uL ,DuL,Dul,DUl,Dul,N  ,N  ,N  ,N  ,N  ,N  ,Dul,DuL,DUL,D  ,INF,INF,
        /*C*/ INF,N  ,N  ,N  ,uL ,uL ,uL ,DuL,Dul,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,DuL,DuL,DuL,DUL,D  ,INF,
        /*C*/ INF,INF,N  ,N  ,uL ,uL ,uL ,DuL,Dul,N  ,N  ,N  ,N  ,N  ,N  ,N  ,N  ,DuL,DuL,DuL,DuL,DUL,D
        }
    };
}();

// Local alignment for proteins (amino acid sequence) with BLOSUM62 score and an extremely wide band.
static auto aa27_01 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        "ALIGATOR"_aa27,
        "GALORA"_aa27,
        align_config
            | align_cfg::scoring{aminoacid_scoring_scheme{aminoacid_similarity_matrix::BLOSUM62}}
            | align_cfg::band{static_band{lower_bound{-10000}, upper_bound{10000}}},
        13,
        "GATOR",
        "GALOR",
        alignment_coordinate{column_index_type{3u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{8u}, row_index_type{5u}},
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

} // namespace seqan3::test::alignment::fixture::local::affine::banded
