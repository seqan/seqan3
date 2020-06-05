// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::search_cfg::detail::result_type.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/type_traits>

#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/search/configuration/detail.hpp>
#include <seqan3/search/search_result.hpp>

namespace seqan3::search_cfg::detail
{

/*!\brief The configuration element storing the configured search result type.
 * \ingroup search_configuration
 *
 * \tparam search_result_t The type of the search result to capture; must be a type specialisation of
 *                         seqan3::search_result.
 *
 * \details
 *
 * Implementation of the search result type configuration element.
 *
 * \see seqan3::search_cfg::result_type
 */
template <typename search_result_t>
//!\cond
    requires seqan3::detail::is_type_specialisation_of_v<search_result_t, search_result>
//!\endcond
struct result_type_tag : public pipeable_config_element<result_type_tag<search_result_t>, seqan3::detail::empty_type>
{
    //!\brief The configured seqan3::search_result type.
    using type = search_result_t;

    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::search_config_id id{seqan3::detail::search_config_id::result_type};
};

/*!\brief Configuration element storing the configured seqan3::search_result for the search algorithm.
 * \ingroup search_configuration
 * \tparam search_result_t The search result type to store; must be a type specialisation of seqan3::search_result.
 *
 * \details
 *
 * This configuration element stores the seqan3::search_result type after configuring the
 * search algorithm with the seqan3::detail::search_configurator.
 * The result type can be accessed via the seqan3::detail::search_traits over the corresponding
 * search configuration type.
 * If the stored search result was not added yet to the search configuration the corresponding
 * result type member will deduce to seqan3::detail::empty_type.
 *
 * \note This configuration element is only added internally during the search configuration and is not intended for
 *       public use.
 */
template <typename search_result_t>
//!\cond
    requires seqan3::detail::is_type_specialisation_of_v<search_result_t, search_result>
//!\endcond
inline constexpr detail::result_type_tag<search_result_t> result_type{};
}  // namespace seqan3::search_cfg::detail
