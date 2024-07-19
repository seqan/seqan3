// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides concepts for the configuration classes.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <concepts>
#include <functional>
#include <type_traits>

#include <seqan3/core/detail/template_inspection.hpp>

namespace seqan3
{
//!\cond
// Forward declarations
struct pipeable_config_element;
//!\endcond
} // namespace seqan3

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// compatibility_table
// ----------------------------------------------------------------------------

/*!\brief Declaration of algorithm specific compatibility table.
 * \ingroup core_configuration
 * \tparam algorithm_id_type The type of the algorithm specific id. Algorithm configurations must maintain
 *                           this table to allow validation checks.
 */
template <typename algorithm_id_type>
inline constexpr std::array<std::array<void *, 0>, 0> compatibility_table{};

// ----------------------------------------------------------------------------
// Concept config_element
// ----------------------------------------------------------------------------

/*!\interface seqan3::detail::config_element <>
 * \brief Concept for an algorithm configuration element.
 * \ingroup core_configuration
 *
 * \extends std::copyable
 * \implements seqan3::pipeable_config_element
 */

/*!\name Requirements for seqan3::detail::config_element
 * \relates seqan3::detail::config_element
 * \brief   You can expect this member on all types that satisfy seqan3::detail::config_element.
 * \{
 */
/*!\var id
 * \brief Algorithm specific static id used for internal validation checks.
 */
//!\}
//!\cond
template <typename config_t>
concept config_element = requires {
    requires std::is_base_of_v<seqan3::pipeable_config_element, config_t>;
    requires std::copyable<config_t>;
    { config_t::id };
};
//!\endcond

/*!\interface seqan3::detail::config_element_pipeable_with <>
 * \brief Concept to check if one configuration element can be combined with another configuration element.
 * \ingroup core_configuration
 *
 * \tparam config1_t The type of the first configuration element.
 * \tparam config2_t The type of the second configuration element.
 *
 * \details
 *
 * This concept is fulfilled if:
 *  * both configurations model seqan3::detail::config_element,
 *  * are defined within the same algorithm configuration domain,
 *  * a seqan3::detail::compatibility_table is defined for the configuration elements and
 *    returns `true` for both configurations.
 */
//!\cond
template <typename config1_t, typename config2_t>
concept config_element_pipeable_with =
    config_element<config1_t> && config_element<config2_t>
    && std::same_as<std::remove_cvref_t<decltype(config1_t::id)>, std::remove_cvref_t<decltype(config2_t::id)>>
    && compatibility_table<std::remove_cvref_t<decltype(config1_t::id)>>[static_cast<int32_t>(config1_t::id)]
                                                                        [static_cast<int32_t>(config2_t::id)];
//!\endcond

} // namespace seqan3::detail

namespace seqan3
{
//!\cond
// Forward declaration.
template <detail::config_element... configs_t>
class configuration;
//!\endcond

/*!\brief Helper variable template to test if a configuration element is combineable with another configuration element
 *        or configuration.
 * \ingroup core_configuration
 *
 * \tparam config1_t Either the type of a configuration element or a configuration.
 * \tparam config2_t Either the type of a configuration element or a configuration.
 *
 * \details
 *
 * This helper variable template checks if `config1_t` fulfills the concept requirements
 * seqan3::detail::config_element_pipeable_with `config2_t`. If `config2_t` is a seqan3::configuration, the check will
 * be expanded to every configuration element contained in the configuration type. Only if `config1_t` is combineable
 * with every element stored inside of the given configuration, this helper variable template evaluates to `true`,
 * otherwise `false`.
 * If `config1_t` is a seqan3::configuration the same applies in combination with the configuration element `config2_t`.
 * If both `config1_t` and `config2_t` are seqan3::configuration types, then the cartesian product between the
 * configuration elements of the first configuration and the second configuration are tested.
 *
 * \noapi{This entity is exposition only!}
 */
template <typename config1_t, typename config2_t>
inline constexpr bool is_config_element_combineable_v = detail::config_element_pipeable_with<config1_t, config2_t>;

//!\cond
// Specialised for config2_t == seqan3::configuration
template <typename config1_t, typename... configs2_t>
inline constexpr bool is_config_element_combineable_v<config1_t, configuration<configs2_t...>> =
    (detail::config_element_pipeable_with<config1_t, configs2_t> && ...);

// Specialised for config1_t == seqan3::configuration
template <typename... configs1_t, typename config2_t>
inline constexpr bool is_config_element_combineable_v<configuration<configs1_t...>, config2_t> =
    (detail::config_element_pipeable_with<configs1_t, config2_t> && ...);

// Specialised for config1_t == seqan3::configuration && config2_t == seqan3::configuration
template <typename... configs1_t, typename... configs2_t>
inline constexpr bool is_config_element_combineable_v<configuration<configs1_t...>, configuration<configs2_t...>> =
    (is_config_element_combineable_v<configs1_t, configuration<configs2_t...>> && ...);
//!\endcond

} // namespace seqan3
