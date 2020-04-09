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
    //!\brief ID for the \ref seqan3::align_cfg::alignment_result_capture "alignment_result_capture" option.
    alignment_result_capture,
    aligned_ends, //!< ID for the \ref seqan3::align_cfg::aligned_ends "aligned_ends" option.
    band,         //!< ID for the \ref seqan3::align_cfg::band "band" option.
    debug,        //!< ID for the \ref seqan3::align_cfg::debug "debug" option.
    gap,          //!< ID for the \ref seqan3::align_cfg::gap "gap" option.
    global,       //!< ID for the \ref seqan3::global_alignment "global alignment" option.
    local,        //!< ID for the \ref seqan3::local_alignment "local alignment" option.
    max_error,    //!< ID for the \ref seqan3::align_cfg::max_error "max_error" option.
    parallel,     //!< ID for the \ref seqan3::align_cfg::parallel "parallel" option.
    result,       //!< ID for the \ref seqan3::align_cfg::result "result" option.
    scoring,      //!< ID for the \ref seqan3::align_cfg::scoring "scoring" option.
    vectorise,    //!< ID for the \ref seqan3::align_cfg::vectorise "vectorise" option.
    SIZE          //!< Represents the number of configuration elements.
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
    {   //0  1  2  3  4  5  6  7  8  9 10 11
        { 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  0: alignment_result_capture
        { 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1}, //  1: aligned_ends
        { 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, //  2: band
        { 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1}, //  3: debug
        { 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1}, //  4: gap
        { 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1}, //  5: global
        { 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1}, //  6: local
        { 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1}, //  7: max_error
        { 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1}, //  8: parallel
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1}, //  9: result
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, // 10: scoring
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0}  // 11: vectorise
    }
};

} // namespace seqan3::detail
