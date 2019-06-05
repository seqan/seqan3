// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
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

namespace seqan3::test::alignment::fixture::global::affine::banded
{

using namespace seqan3;
using namespace seqan3::detail;

inline constexpr auto align_config = align_cfg::mode{global_alignment} |
                                     align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
                                     align_cfg::band{static_band{lower_bound{-3}, upper_bound{8}}};

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
        align_config | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        -18,
        "A---ACCGGTTAACCGGTT",
        "ACGTAC----------GTA",
        alignment_coordinate{detail::column_index_type{0u}, detail::row_index_type{0u}},
        alignment_coordinate{detail::column_index_type{16u}, detail::row_index_type{9u}}
    };
}();

} // namespace seqan3::test::alignment::fixture::global::affine::banded
