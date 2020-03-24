// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <vector>

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_mode.hpp>
#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "alignment_fixture.hpp"

using seqan3::operator""_dna4;

namespace seqan3::test::alignment::fixture::global::affine::banded
{

inline constexpr auto align_config = seqan3::align_cfg::mode{seqan3::global_alignment} |
                                     seqan3::align_cfg::gap{seqan3::gap_scheme{seqan3::gap_score{-1},
                                                                               seqan3::gap_open_score{-10}}} |
                                     seqan3::align_cfg::band{seqan3::static_band{seqan3::lower_bound{-3},
                                                                                 seqan3::upper_bound{8}}};

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
        align_config | seqan3::align_cfg::scoring{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4},
                                                                                    seqan3::mismatch_score{-5}}},
        -18,
        "A---ACCGGTTAACCGGTT",
        "ACGTAC----------GTA",
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{0u}, seqan3::detail::row_index_type{0u}},
        seqan3::alignment_coordinate{seqan3::detail::column_index_type{16u}, seqan3::detail::row_index_type{9u}},
        std::vector<std::optional<int32_t>>
        {
        //       A  ,A  ,C  ,C  ,G  ,G  ,T  ,T  ,A  ,A  ,C  ,C  ,G  ,G  ,T  ,T  ,
             0  ,-11,-12,-13,-14,-15,-16,-17,-18,INF,INF,INF,INF,INF,INF,INF,INF,
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
        //       A  ,A  ,C  ,C  ,G  ,G  ,T  ,T  ,A  ,A  ,C  ,C  ,G  ,G  ,T  ,T  ,
             N  ,L  ,l  ,l  ,l  ,l  ,l  ,l  ,l  ,INF,INF,INF,INF,INF,INF,INF,INF,
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

} // namespace seqan3::test::alignment::fixture::global::affine::banded
