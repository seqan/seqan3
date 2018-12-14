// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides compatibility matrix for search configurations.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/configuration_utility.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// search_config_id
// ----------------------------------------------------------------------------

/*!\brief Specifies an id for every configuration element.
 * \ingroup configuration
 *
 * \details
 *
 * The seqan3::search_cfg::id is used to identify a specific search configuration element independent of
 * it's concrete type and position within the seqan3::search_cfg::search_configuration object.
 * Thus one can access the value of the corresponding configuration element via the special get interface.
 *
 * ### Example
 *
 * ```cpp
 * search_cfg::search_configuration cfg = search_cfg::max_total_errors(3);
 * auto max_total_errors = get<search_cfg::id::max_total_errors>(cfg);  // max_total_errors = 3;
 * ```
 */
enum struct search_config_id : uint8_t
{
    //!\brief Identifier for max_errors configuration.
    max_error,
    max_error_rate,
    output,
    mode,
    //!\cond
    // ATTENTION: Must always be the last item; will be used to determine the number of ids.
    SIZE
    //!\endcond
};

// ----------------------------------------------------------------------------
// search_config_validation_matrix
// ----------------------------------------------------------------------------

/*!\brief Compatibility matrix to check how search configuration elements can be combined.
 * \ingroup configuration
 *
 * \details
 *
 * This matrix is used to check if the specified search configurations can be combined with each other.
 * A cell value `true`, indicates that the corresponding seqan3::detail::search_config_id in the current column can
 * be combined with the associated seqan3::search_cfg::id in the current row. The size of the matrix is determined by
 * the enum value `SIZE` of seqan3::detail::search_config_id.
 */
template <>
inline constexpr std::array<std::array<bool, static_cast<uint8_t>(search_config_id::SIZE)>,
                            static_cast<uint8_t>(search_config_id::SIZE)> compatibility_table<search_config_id> =
{
    {
        // max_error, max_error_rate, output, mode
        { 0, 0, 1, 1 },
        { 0, 0, 1, 1 },
        { 1, 1, 0, 1 },
        { 1, 1, 1, 0 }
    }
};

} // namespace seqan3::detail
