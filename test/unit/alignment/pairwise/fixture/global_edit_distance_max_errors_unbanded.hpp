
#pragma once

#include <vector>

#include "alignment_fixture.hpp"
#include "global_edit_distance_unbanded.hpp"

#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_max_error.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

/**
 * NOTE: max_errors is a special case where it will produces the same matrix excepts that it will cutoff all scores from
 * the bottom to the top in the matrix until the score does not exceed the allowed error anymore.
 *
 * Thus we can apply some masking to the matrix and get similiar
 */

namespace seqan3::test::alignment::fixture::global::edit_distance::max_errors::unbanded
{

using namespace seqan3::test::alignment::fixture::global::edit_distance::unbanded;
using detail::column_index_type;
using detail::row_index_type;

static auto dna4_01_e255 = []()
{
    return alignment_fixture
    {
        "AACCGGTTAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        align_cfg::edit | align_cfg::max_error{255},
        -8,
        "AACCGGTTAACCGGTT",
        "A-C-G-T-A-C-G-TA",
        dna4_01.front_coordinate,
        dna4_01.back_coordinate,
        dna4_01.score_vector,
        dna4_01.trace_vector
    };
}();

static auto dna4_01T_e255 = []()
{
    return alignment_fixture
    {
        "ACGTACGTA"_dna4,
        "AACCGGTTAACCGGTT"_dna4,
        align_cfg::edit | align_cfg::max_error{255},
        -8,
        "A-C-G-T-A-C-G-TA",
        "AACCGGTTAACCGGTT",
        dna4_01T.front_coordinate,
        dna4_01T.back_coordinate,
        dna4_01T.score_vector,
        dna4_01T.trace_vector
    };
}();

static auto dna4_02_e255 = []()
{
    return alignment_fixture
    {
        "AACCGGTAAACCGGTT"_dna4,
        "ACGTACGTA"_dna4,
        align_cfg::edit | align_cfg::max_error{255},
        -8,
        "AACCGGTAAACCGGTT",
        "A-C-G-TA--C-G-TA",
        dna4_02.front_coordinate,
        dna4_02.back_coordinate,
        dna4_02.score_vector,
        dna4_02.trace_vector
    };
}();

static auto aa27_01_e255 = []()
{
    return alignment_fixture
    {
        "UUWWRRIIUUWWRRII"_aa27,
        "UWRIUWRIU"_aa27,
        align_cfg::edit | align_cfg::max_error{255},
        -8,
        "UUWWRRIIUUWWRRII",
        "U-W-R-I-U-W-R-IU",
        aa27_01.front_coordinate,
        aa27_01.back_coordinate,
        aa27_01.score_vector,
        aa27_01.trace_vector
    };
}();

static auto aa27_01T_e255 = []()
{
    return alignment_fixture
    {
        "UWRIUWRIU"_aa27,
        "UUWWRRIIUUWWRRII"_aa27,
        align_cfg::edit | align_cfg::max_error{255},
        -8,
        "U-W-R-I-U-W-R-IU",
        "UUWWRRIIUUWWRRII",
        aa27_01T.front_coordinate,
        aa27_01T.back_coordinate,
        aa27_01T.score_vector,
        aa27_01T.trace_vector
    };
}();

} // namespace seqan3::test::alignment::fixture::global::edit_distance::unbanded
