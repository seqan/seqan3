// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <vector>

#include <seqan3/alignment/configuration/align_config_mode.hpp>
#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>

#include "alignment_fixture.hpp"

namespace seqan3::test::alignment::fixture::local::affine::unbanded
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
        align_config | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        11,
        "GTTTA",
        "GTCTA",
        alignment_coordinate{column_index_type{5u}, row_index_type{2u}},
        alignment_coordinate{column_index_type{10u}, row_index_type{7u}}
    };
}();

// The same alignment with sequences swapped.
static auto dna4_02 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        "ACGTCTACGTA"_dna4,
        "AACCGGTTTAACCGGTT"_dna4,
        align_config | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        11,
        "GTCTA",
        "GTTTA",
        alignment_coordinate{column_index_type{2u}, row_index_type{5u}},
        alignment_coordinate{column_index_type{7u}, row_index_type{10u}}
    };
}();

// Local alignment starting in the first row. Verifies that free end gaps are performed correctly.
static auto dna4_03 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        "ataagcgtctcg"_dna4,
        "tcatagagttgc"_dna4,
        align_cfg::mode{local_alignment}
            | align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-1}}}
            | align_cfg::scoring{nucleotide_scoring_scheme{match_score{2}, mismatch_score{-1}}},
        9,
        "ATAAGCGT",
        "AT-AGAGT",
        alignment_coordinate{column_index_type{0u}, row_index_type{2u}},
        alignment_coordinate{column_index_type{8u}, row_index_type{9u}}
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
        align_config | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        0,
        "",
        "",
        alignment_coordinate{column_index_type{0u}, row_index_type{6u}}, // in SeqAn2 both coordinates are (0,0)
        alignment_coordinate{column_index_type{0u}, row_index_type{6u}}
    };
}();

// Local alignment in the begin and end of sequences.
static auto dna4_05 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        "AAAAAATCCCCCC"_dna4,
        "CCCCCCTAAAAAA"_dna4,
        align_config | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        24,
        "AAAAAA",
        "AAAAAA",
        alignment_coordinate{column_index_type{0u}, row_index_type{7u}},
        alignment_coordinate{column_index_type{6u}, row_index_type{13u}}
    };
}();

// Local RNA alignment with a longer sequence of gaps.
static auto rna5_01 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        "AAAAAAUUUUNNUUUUCCCCCC"_rna5,
        "AAAAAACCCCCC"_rna5,
        align_config | align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}},
        28,
        "AAAAAAUUUUNNUUUUCCCCCC",
        "AAAAAA----------CCCCCC",
        alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{22u}, row_index_type{12u}}
    };
}();

// Local alignment for proteins (amino acid sequence) with BLOSUM62 score.
static auto aa27_01 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        "ALIGATOR"_aa27,
        "GALORA"_aa27,
        align_config | align_cfg::scoring{aminoacid_scoring_scheme{aminoacid_similarity_matrix::BLOSUM62}},
        13,
        "GATOR",
        "GALOR",
        alignment_coordinate{column_index_type{3u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{8u}, row_index_type{5u}}
    };
}();

// Local alignment with empty sequence.
static auto aa27_02 = []()
{
    using detail::column_index_type;
    using detail::row_index_type;

    return alignment_fixture
    {
        "ALIGATOR"_aa27,
        ""_aa27,
        align_config | align_cfg::scoring{aminoacid_scoring_scheme{aminoacid_similarity_matrix::BLOSUM62}},
        0,
        "",
        "",
        alignment_coordinate{column_index_type{0u}, row_index_type{0u}},
        alignment_coordinate{column_index_type{0u}, row_index_type{0u}}
    };
}();

} // namespace seqan3::test::alignment::fixture::local::affine::unbanded
