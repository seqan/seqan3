// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::search_traits.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/hit.hpp>
#include <seqan3/search/configuration/output.hpp>
#include <seqan3/search/configuration/result_type.hpp>

namespace seqan3::detail
{

/*!\brief A collection of traits extracted from the search configuration.
 * \ingroup search
 *
 * \tparam search_configuration_t The type of the search algorithm configuration; must be of type seqan3::configuration.
 */
template <typename search_configuration_t>
struct search_traits
{
    //!\brief A specific empty search result type indicating misconfiguration of the search algorithm.
    using empty_search_result_type = search_result<empty_type, empty_type, empty_type, empty_type>;
    //!\brief The configured search result type.
    using search_result_type =
        typename remove_cvref_t<decltype(std::declval<search_configuration_t>().get_or(
                search_cfg::detail::result_type<empty_search_result_type>))
        >::type;

    //!\brief A flag indicating whether search should be invoked with total errors.
    static constexpr bool has_max_error_total =
                                        search_configuration_t::template exists<search_cfg::max_error_total>();
    //!\brief A flag indicating whether search should be invoked with substitution errors.
    static constexpr bool has_max_error_substitution =
                                        search_configuration_t::template exists<search_cfg::max_error_substitution>();
    //!\brief A flag indicating whether search should be invoked with insertion errors.
    static constexpr bool has_max_error_insertion =
                                        search_configuration_t::template exists<search_cfg::max_error_insertion>();
    //!\brief A flag indicating whether search should be invoked with deletion errors.
    static constexpr bool has_max_error_deletion =
                                        search_configuration_t::template exists<search_cfg::max_error_deletion>();

    //!\brief A flag that indicates whether the search should be invoked with only specified total errors.
    static constexpr bool only_max_error_total = has_max_error_total &&
                                                 !has_max_error_substitution &&
                                                 !has_max_error_insertion &&
                                                 !has_max_error_deletion;

    //!\brief A flag indicating whether search should find all hits.
    static constexpr bool search_all_hits = search_configuration_t::template exists<hit_all_tag>();
    //!\brief A flag indicating whether search should find best hits.
    static constexpr bool search_single_best_hit = search_configuration_t::template exists<hit_single_best_tag>();
    //!\brief A flag indicating whether search should find all best hits.
    static constexpr bool search_all_best_hits = search_configuration_t::template exists<hit_all_best_tag>();
    //!\brief A flag indicating whether search should find strata hits.
    static constexpr bool search_strata_hits = search_configuration_t::template exists<search_cfg::hit_strata>();
    //!\brief A flag indicating whether hit configuration was set in the search configuration.
    static constexpr bool has_hit_configuration = search_all_hits ||
                                                  search_single_best_hit ||
                                                  search_all_best_hits ||
                                                  search_strata_hits ||
                                                  search_configuration_t::template exists<search_cfg::hit>();

    //!\brief A flag indicating whether search should return the query_id.
    static constexpr bool output_query_id = search_configuration_t::template exists<detail::output_query_id_tag>();
    //!\brief A flag indicating whether search should return the reference_id.
    static constexpr bool output_reference_id =
                              search_configuration_t::template exists<detail::output_reference_id_tag>();
    //!\brief A flag indicating whether search should return the reference_begin_position.
    static constexpr bool output_reference_begin_position =
                              search_configuration_t::template exists<detail::output_reference_begin_position_tag>();
    //!\brief A flag indicating whether search should return the index_cursor.
    static constexpr bool output_index_cursor =
                              search_configuration_t::template exists<detail::output_index_cursor_tag>();
    //!\brief A flag indicating whether it is required to call cursor.locate() to retrieve the respective information.
    static constexpr bool output_requires_locate_call = output_reference_id | output_reference_begin_position;

    //!\brief A flag indicating whether output configuration was set in the search configuration.
    static constexpr bool has_output_configuration = output_query_id |
                                                     output_reference_id |
                                                     output_reference_begin_position |
                                                     output_index_cursor;
};

} // namespace seqan3::detail
