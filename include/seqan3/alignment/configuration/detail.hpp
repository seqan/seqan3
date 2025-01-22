// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides some utility functions for the alignment configurations.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/configuration/detail/concept.hpp>

namespace seqan3::detail
{

/*!\brief An internal enum to check for a consistent configuration object.
 * \ingroup alignment_configuration
 */
enum struct align_config_id : uint8_t
{
    band,                  //!< ID for the \ref seqan3::align_cfg::band_fixed_size "band" option.
    debug,                 //!< ID for the \ref seqan3::align_cfg::detail::debug "debug" option.
    gap,                   //!< ID for the \ref seqan3::align_cfg::gap_cost_affine "gap_cost_affine" option.
    global,                //!< ID for the \ref seqan3::align_cfg::method_global "global alignment" option.
    local,                 //!< ID for the \ref seqan3::align_cfg::method_local "local alignment" option.
    min_score,             //!< ID for the \ref seqan3::align_cfg::min_score "min_score" option.
    on_result,             //!< ID for the \ref seqan3::align_cfg::on_result "on_result" option.
    output_alignment,      //!< ID for the \ref seqan3::align_cfg::output_alignment "alignment output" option.
    output_begin_position, //!< ID for the \ref seqan3::align_cfg::output_begin_position "begin position output" option.
    output_end_position,   //!< ID for the \ref seqan3::align_cfg::output_end_position "end position output" option.
    output_sequence1_id,   //!< ID for the \ref seqan3::align_cfg::output_sequence1_id "sequence1 id output" option.
    output_sequence2_id,   //!< ID for the \ref seqan3::align_cfg::output_sequence2_id "sequence2 id output" option.
    output_score,          //!< ID for the \ref seqan3::align_cfg::output_score "score output" option.
    parallel,              //!< ID for the \ref seqan3::align_cfg::parallel "parallel" option.
    result_type,           //!< ID for the \ref seqan3::align_cfg::detail::result_type "result_type" option.
    score_type,            //!< ID for the \ref seqan3::align_cfg::score_type "score_type" option.
    scoring,               //!< ID for the \ref seqan3::align_cfg::scoring_scheme "scoring_scheme" option.
    vectorised,            //!< ID for the \ref seqan3::align_cfg::vectorised "vectorised" option.
    SIZE                   //!< Represents the number of configuration elements.
};

// ----------------------------------------------------------------------------
// compatibility_table
// ----------------------------------------------------------------------------

/*!\brief Declaration of algorithm specific compatibility table.
 * \ingroup alignment_configuration
 * \tparam algorithm_id_type The type of the algorithm specific id. Algorithm configurations must maintain
 *                           this table to allow validation checks.
 */
template <>
inline constexpr std::array<std::array<bool, static_cast<uint8_t>(align_config_id::SIZE)>,
                            static_cast<uint8_t>(align_config_id::SIZE)>
    compatibility_table<align_config_id>{{
        //band
        //|  debug
        //|  |  gap
        //|  |  |  global
        //|  |  |  |  local
        //|  |  |  |  |  min_score
        //|  |  |  |  |  |  on_result
        //|  |  |  |  |  |  |  output_alignment
        //|  |  |  |  |  |  |  |  output_begin_position
        //|  |  |  |  |  |  |  |  |  output_end_position
        //|  |  |  |  |  |  |  |  |  |  output_sequence1_id
        //|  |  |  |  |  |  |  |  |  |  |  output_sequence2_id
        //|  |  |  |  |  |  |  |  |  |  |  |  output_score
        //|  |  |  |  |  |  |  |  |  |  |  |  |  parallel
        //|  |  |  |  |  |  |  |  |  |  |  |  |  |  result_type
        //|  |  |  |  |  |  |  |  |  |  |  |  |  |  | score_type
        //|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  scoring
        //|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  vectorised
        {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  0: band
        {1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  1: debug
        {1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  2: gap
        {1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  3: global
        {1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  4: local
        {1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  5: max_error
        {1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  6: on_result
        {1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  7: output_alignment
        {1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  8: output_begin_position
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1}, // 9: output_end_position
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1}, // 10: output_sequence1_id
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1}, // 11: output_sequence2_id
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1}, // 12: output_score
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1}, // 13: parallel
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1}, // 14: result_type
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1}, // 15: score_type
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, // 16: scoring
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0}  // 17: vectorised
    }};

} // namespace seqan3::detail
