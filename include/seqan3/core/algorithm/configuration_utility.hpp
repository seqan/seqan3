// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides functionality to access get function by enum values.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <type_traits>

#include <meta/meta.hpp>

#include <seqan3/core/algorithm/concept.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/core/concept/tuple.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// compatibility_table
// ----------------------------------------------------------------------------

/*!\brief Declaration of algorithm specific compatibility table.
 * \ingroup algorithm
 * \tparam algorithm_id_type The type of the algorithm specific id. Algorithm configurations must maintain
 *                           this table to allow validation checks.
 */
template <typename algorithm_id_type>
inline constexpr std::array<std::array<void, 0>, 0> compatibility_table;

// ----------------------------------------------------------------------------
// Metafunction is_configuration_valid
// ----------------------------------------------------------------------------

/*!\brief Value metafunction which checks if a given type is compatible with a list of other types.
 * \ingroup algorithm
 * \tparam query_t       The type to check for compatibility.
 * \tparam compare_types The types to compare against.
 *
 * \details
 *
 * Checks if the type is from the same algorithm configuration and if it can be combined with any of the
 * existing elements in the current configuration.
 *
 * \see seqan3::detail::is_configuration_valid_v
 */
template <config_element_concept query_t, config_element_concept ... compare_types>
struct is_configuration_valid :
    public std::conditional_t<
        (std::is_same_v<remove_cvref_t<decltype(query_t::id)>, remove_cvref_t<decltype(compare_types::id)>> && ...) &&
        (compatibility_table<remove_cvref_t<decltype(query_t::id)>>
                [static_cast<std::underlying_type_t<remove_cvref_t<decltype(query_t::id)>>>(query_t::id)]
                [static_cast<std::underlying_type_t<remove_cvref_t<decltype(query_t::id)>>>(compare_types::id)] && ...),
        std::true_type,  // If condition is true.
        std::false_type  // If condition is false.
    >
{};

/*!\brief Helper variable template to check for valid configuration compositions.
 * \ingroup algorithm
 * \see seqan3::detail::is_configuration_valid
 */
template <typename query_t, typename ... compare_types>
inline constexpr bool is_configuration_valid_v = is_configuration_valid<query_t, compare_types...>::value;

// ----------------------------------------------------------------------------
// Metafunction is_same_configuration_f
// ----------------------------------------------------------------------------

/*!\brief Helper meta function to check if a template type is contained in a seqan3::configuration.
 * \ingroup algorithm
 *
 * \details
 *
 * This helper meta function is used to provide the `get` and `value_or` interface for template template types.
 */
template <template <typename ...> typename query_t>
struct is_same_configuration_f
{
    /*!\brief A type template that evaluates to std::true_type if the given type is a specialization of `query_t`,
     *        otherwise std::false_type.
     * \tparam compare_type The type to compare against `query_t`.
     */
    template <typename compare_type>
    using invoke = is_type_specialisation_of<compare_type, query_t>;
};

// ----------------------------------------------------------------------------
// Metafunction is_algorithm_configuration
// ----------------------------------------------------------------------------

/*!\brief Value metafunction that returns whether a type is an algorithm configuration.
 * \ingroup algorithm
 *
 * \returns std::true_type if the given type is a seqan3::detail::configuration, else std::false_type.
 */
template <typename object_t>
struct is_algorithm_configuration : std::false_type
{};

//!\cond
template <typename ...config_elements_t>
struct is_algorithm_configuration<seqan3::configuration<config_elements_t...>> : std::true_type
{};
//!\endcond

/*!\brief Helper variable template for seqan3::detail::is_algorithm_configuration.
 * \ingroup algorithm
 */
template <typename object_t>
inline constexpr bool is_algorithm_configuration_v = is_algorithm_configuration<object_t>::value;

} // namespace seqan3::detail
