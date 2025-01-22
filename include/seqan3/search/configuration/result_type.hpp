// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::search_cfg::detail::result_type.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/configuration/pipeable_config_element.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/search/configuration/detail.hpp>
#include <seqan3/search/search_result.hpp>

namespace seqan3::search_cfg::detail
{

/*!\brief Configuration element storing the configured seqan3::search_result for the search algorithm.
 * \ingroup search_configuration
 * \see search_configuration
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
    requires seqan3::detail::is_type_specialisation_of_v<search_result_t, search_result>
class result_type : private pipeable_config_element
{
public:
    //!\brief The configured seqan3::search_result type.
    using type = search_result_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr result_type() = default;                                //!< Defaulted.
    constexpr result_type(result_type const &) = default;             //!< Defaulted.
    constexpr result_type(result_type &&) = default;                  //!< Defaulted.
    constexpr result_type & operator=(result_type const &) = default; //!< Defaulted.
    constexpr result_type & operator=(result_type &&) = default;      //!< Defaulted.
    ~result_type() = default;                                         //!< Defaulted.

    //!\}

    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::search_config_id id{seqan3::detail::search_config_id::result_type};
};
} // namespace seqan3::search_cfg::detail
