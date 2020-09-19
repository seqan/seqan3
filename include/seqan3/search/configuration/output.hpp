// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the configuration for the content of the search result.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/search/configuration/detail.hpp>

namespace seqan3::search_cfg
{

/*!\brief Include the query_id in the seqan3::search_result returned by a call to seqan3::search.
 * \ingroup search_configuration
 * \sa \ref search_configuration_subsection_output "Section on Output"
 */
struct output_query_id : public pipeable_config_element<output_query_id>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::output_query_id};
};

/*!\brief Include the reference_id in the seqan3::search_result returned by a call to seqan3::search.
 * \ingroup search_configuration
 * \sa \ref search_configuration_subsection_output "Section on Output"
 */
struct output_reference_id : public pipeable_config_element<output_reference_id>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::output_reference_id};
};

/*!\brief Include the reference_begin_position in the seqan3::search_result returned by a call to seqan3::search.
 * \ingroup search_configuration
 * \sa \ref search_configuration_subsection_output "Section on Output"
 */
struct output_reference_begin_position : public pipeable_config_element<output_reference_begin_position>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::output_reference_begin_position};
};

/*!\brief Include the index_cursor in the seqan3::search_result returned by a call to seqan3::search.
 * \ingroup search_configuration
 * \sa \ref search_configuration_subsection_output "Section on Output"
 */
struct output_index_cursor : public pipeable_config_element<output_index_cursor>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::output_index_cursor};
};

} // namespace seqan3::search_cfg
