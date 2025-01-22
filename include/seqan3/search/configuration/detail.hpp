// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides compatibility matrix for search configurations.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/configuration/detail/concept.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// search_config_id
// ----------------------------------------------------------------------------

/*!\brief Specifies an id for every configuration element.
 * \ingroup search_configuration
 * \see search_configuration
 *
 * \details
 *
 * The seqan3::detail::search_config_id used to identify a specific search configuration element independent of
 * its concrete type and position within the \ref seqan3::search_cfg "search configuration object".
 * Thus one can access the value of the corresponding configuration element via the special get interface.
 *
 * ### Example
 *
 * ```cpp
 * search_cfg cfg = search_cfg::max_total_errors(3);
 * auto max_total_errors = get<search_cfg::id::max_total_errors>(cfg);  // max_total_errors = 3;
 * ```
 */
enum struct search_config_id : uint8_t
{
    max_error_total,                 //!< Identifier for the max_error_total configuration.
    max_error_substitution,          //!< Identifier for the max_error_substitution configuration.
    max_error_insertion,             //!< Identifier for the max_error_insertion configuration.
    max_error_deletion,              //!< Identifier for the max_error_deletion configuration.
    on_result,                       //!< Identifier for the configuration to pass a user defined callable.
    output_query_id,                 //!< Identifier for the output configuration of the query_id.
    output_reference_id,             //!< Identifier for the output configuration of the reference_id.
    output_reference_begin_position, //!< Identifier for the output configuration of the reference_begin_position.
    output_index_cursor,             //!< Identifier for the output configuration of the index_cursor.
    hit,                             //!< Identifier for the hit configuration (all, all_best, single_best, strata).
    parallel,                        //!< Identifier for the parallel execution configuration.
    result_type,                     //!< Identifier for the configured search result type.
    //!\cond
    // ATTENTION: Must always be the last item; will be used to determine the number of ids.
    SIZE //!< Determines the size of the enum.
    //!\endcond
};

// ----------------------------------------------------------------------------
// search_config_validation_matrix
// ----------------------------------------------------------------------------

/*!\brief Compatibility matrix to check how search configuration elements can be combined.
 * \ingroup search_configuration
 * \see search_configuration
 *
 * \details
 *
 * This matrix is used to check if the specified search configurations can be combined with each other.
 * A cell value `true`, indicates that the corresponding seqan3::detail::search_config_id in the current column can
 * be combined with the associated seqan3::detail::search_config in the current row. The size of the matrix is
 * determined by the enum value `SIZE` of seqan3::detail::search_config_id.
 */
template <>
inline constexpr std::array<std::array<bool, static_cast<uint8_t>(search_config_id::SIZE)>,
                            static_cast<uint8_t>(search_config_id::SIZE)>
    compatibility_table<search_config_id> = {{
        // max_error_total,
        // |  max_error_substitution,
        // |  |  max_error_insertion,
        // |  |  |  max_error_deletion,
        // |  |  |  |  on_result,
        // |  |  |  |  |  output_query_id,
        // |  |  |  |  |  |  output_reference_id,
        // |  |  |  |  |  |  |  output_reference_begin_position,
        // |  |  |  |  |  |  |  |  output_index_cursor,
        // |  |  |  |  |  |  |  |  |  hit,
        // |  |  |  |  |  |  |  |  |  |  parallel,
        // |  |  |  |  |  |  |  |  |  |  |  result_type
        {0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, // max_error_total
        {1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, // max_error_substitution
        {1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1}, // max_error_insertion
        {1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1}, // max_error_deletion
        {1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1}, // on_result
        {1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1}, // output_query_id
        {1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1}, // output_reference_id
        {1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1}, // output_reference_begin_position
        {1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1}, // output_index_cursor
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1}, // hit
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1}, // parallel
        {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0}  // result_type
    }};

} // namespace seqan3::detail
