// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides some utility functions for the alignment configurations.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/configuration_utility.hpp>

namespace seqan3::detail
{

/*!\brief An internal enum to check for a consistent configuration object.
 * \ingroup alignment_configuration
 */
enum struct align_config_id : uint8_t
{

    aligned_ends,          //!< ID for the \ref seqan3::align_cfg::aligned_ends "aligned_ends" option.
    band,                  //!< ID for the \ref seqan3::align_cfg::band_fixed_size "band" option.
    debug,                 //!< ID for the \ref seqan3::align_cfg::detail::debug "debug" option.
    gap,                   //!< ID for the \ref seqan3::align_cfg::gap "gap" option.
    global,                //!< ID for the \ref seqan3::align_cfg::method_global "global alignment" option.
    local,                 //!< ID for the \ref seqan3::align_cfg::method_local "local alignment" option.
    max_error,             //!< ID for the \ref seqan3::align_cfg::max_error "max_error" option.
    on_result,             //!< ID for the \ref seqan3::align_cfg::on_result "on_result" option.
    output_alignment,      //!< ID for the \ref seqan3::align_cfg::output_alignment "alignment output" option.
    output_begin_position, //!< ID for the \ref seqan3::align_cfg::output_begin_position "begin position output" option.
    output_end_position,   //!< ID for the \ref seqan3::align_cfg::output_end_position "end position output" option.
    output_sequence1_id,   //!< ID for the \ref seqan3::align_cfg::output_sequence1_id "sequence1 id output" option.
    output_sequence2_id,   //!< ID for the \ref seqan3::align_cfg::output_sequence2_id "sequence2 id output" option.
    output_score,          //!< ID for the \ref seqan3::align_cfg::output_score "score output" option.
    parallel,              //!< ID for the \ref seqan3::align_cfg::parallel "parallel" option.
    result,                //!< ID for the \ref seqan3::align_cfg::result "result" option.
    result_type,           //!< ID for the \ref seqan3::align_cfg::detail::result_type "result_type" option.
    scoring,               //!< ID for the \ref seqan3::align_cfg::scoring_scheme "scoring_scheme" option.
    vectorised,            //!< ID for the \ref seqan3::align_cfg::vectorised "vectorised" option.
    SIZE                   //!< Represents the number of configuration elements.
};

// ----------------------------------------------------------------------------
// compatibility_table
// ----------------------------------------------------------------------------

/*!\brief Declaration of algorithm specific compatibility table.
 * \ingroup algorithm
 * \tparam algorithm_id_type The type of the algorithm specific id. Algorithm configurations must maintain
 *                           this table to allow validation checks.
 */
template <>
inline constexpr std::array<std::array<bool, static_cast<uint8_t>(align_config_id::SIZE)>,
                            static_cast<uint8_t>(align_config_id::SIZE)> compatibility_table<align_config_id>
{
    {   //aligned_ends
        //|  band
        //|  |  debug
        //|  |  |  gap
        //|  |  |  |  global
        //|  |  |  |  |  local
        //|  |  |  |  |  |  max_error
        //|  |  |  |  |  |  |  on_result
        //|  |  |  |  |  |  |  |  output_alignment
        //|  |  |  |  |  |  |  |  |  output_begin_position
        //|  |  |  |  |  |  |  |  |  |  output_end_position
        //|  |  |  |  |  |  |  |  |  |  |  output_sequence1_id
        //|  |  |  |  |  |  |  |  |  |  |  |  output_sequence2_id
        //|  |  |  |  |  |  |  |  |  |  |  |  |  output_score
        //|  |  |  |  |  |  |  |  |  |  |  |  |  |  parallel
        //|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  result
        //|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  result_type
        //|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  scoring
        //|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  vectorised
        { 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  0: aligned_ends
        { 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  1: band
        { 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  2: debug
        { 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  3: gap
        { 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  4: global
        { 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  5: local
        { 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  6: max_error
        { 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  7: on_result
        { 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  8: output_alignment
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  9: output_begin_position
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1}, // 10: output_end_position
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1}, // 11: output_sequence1_id
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1}, // 12: output_sequence2_id
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1}, // 13: output_score
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1}, // 14: parallel
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1}, // 15: result
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1}, // 16: result_type
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, // 17: scoring
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0}  // 18: vectorised
    }
};

} // namespace seqan3::detail
