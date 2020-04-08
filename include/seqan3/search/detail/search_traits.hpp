// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::search_traits.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/max_error_rate.hpp>
#include <seqan3/search/configuration/mode.hpp>
#include <seqan3/search/configuration/output.hpp>

namespace seqan3::detail
{

/*!\brief A collection of traits extracted from the search configuration.
 * \ingroup search
 *
 * \tparam search_configuration_t The type of the search algorithm configuration; must be of type
 *                                seqan3::configuration.
 */
template <typename search_configuration_t>
//!\cond
    requires is_type_specialisation_of_v<remove_cvref_t<search_configuration_t>, configuration>
//!\endcond
struct search_traits
{
    //!\brief A flag indicating whether search should be invoked with max error.
    static constexpr bool search_with_max_error = search_configuration_t::template exists<search_cfg::max_error>();
    //!\brief A flag indicating whether search should be invoked with max error rate.
    static constexpr bool search_with_max_error_rate =
        search_configuration_t::template exists<search_cfg::max_error_rate>();
    //!\brief A flag indicating whether error configuration was set in the search configuration.
    static constexpr bool has_error_configuration = search_with_max_error | search_with_max_error_rate;

    //!\brief A flag indicating whether search should find all hits.
    static constexpr bool search_all_hits =
        search_configuration_t::template exists<search_cfg::mode<detail::search_mode_all>>();
    //!\brief A flag indicating whether search should find best hits.
    static constexpr bool search_best_hits =
        search_configuration_t::template exists<search_cfg::mode<detail::search_mode_best>>();
    //!\brief A flag indicating whether search should find all best hits.
    static constexpr bool search_all_best_hits =
        search_configuration_t::template exists<search_cfg::mode<detail::search_mode_all_best>>();
    //!\brief A flag indicating whether search should find strata hits.
    static constexpr bool search_strata_hits =
        search_configuration_t::template exists<search_cfg::mode<search_cfg::strata>>();
    //!\brief A flag indicating whether mode configuration was set in the search configuration.
    static constexpr bool has_mode_configuration = search_all_hits |
                                                   search_best_hits |
                                                   search_all_best_hits |
                                                   search_strata_hits;

    //!\brief A flag indicating whether search should return the index cursor.
    static constexpr bool search_return_index_cursor =
        search_configuration_t::template exists<search_cfg::output<detail::search_output_index_cursor>>();
    //!\brief A flag indicating whether search should return the text position.
    static constexpr bool search_return_text_position =
        search_configuration_t::template exists<search_cfg::output<detail::search_output_text_position>>();
    //!\brief A flag indicating whether output configuration was set in the search configuration.
    static constexpr bool has_output_configuration = search_return_index_cursor | search_return_text_position;
};

} // namespace seqan3::detail
