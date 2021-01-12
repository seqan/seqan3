// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concepts for the configuration classes.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concepts>
#include <functional>
#include <type_traits>

#include <meta/meta.hpp>

#include <seqan3/core/detail/template_inspection.hpp>

namespace seqan3
{
//!\cond
// Forward declarations
template <typename derived_t, typename value_t = void>
struct pipeable_config_element;
//!\endcond
}

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
inline constexpr std::array<std::array<void *, 0>, 0> compatibility_table{};

// ----------------------------------------------------------------------------
// Concept config_element
// ----------------------------------------------------------------------------

#if SEQAN3_WORKAROUND_GCC_PIPEABLE_CONFIG_CONCEPT
/*!\brief A helper class to check if a type has a static member called `id`.
 * \ingroup algorithm
 *
 * \details
 *
 * This class is needed for gcc versions older than 11. It adds a SFINAE check to test if a type
 * has a static member called `id`, which is needed for the concept defintions around the
 * pipeable configuration element concepts.
 */
struct config_id_accessor
{
private:
    //!\brief Helper variable template to convert the configuration enum identifier to an integer.
    template <typename config_t>
    static constexpr int32_t as_int = static_cast<int32_t>(std::remove_cvref_t<config_t>::id);

    //!\brief Helper function to check if static id member exists.
    template <typename config_t>
    static constexpr auto has_id_member(int) -> decltype((static_cast<void>(config_t::id), true))
    {
        return true;
    }

    //!\overload
    template <typename config_t>
    static constexpr bool has_id_member(...)
    {
        return false;
    }

public:
    //!\brief Type alias for the internal enumeration type that represents the id.
    template <typename config_t>
    using id_type = std::remove_cvref_t<decltype(config_t::id)>;

    /*!\brief Checks if two configuration types are compatible.
     *
     * \tparam config1_t The type of the first configuration to check against the second.
     * \tparam config2_t The type of the second configuration.
     *
     * \details
     *
     * Uses the seqan3::detail::compatibility_table of the corresponding configuration element domain to check
     * if both configurations can be combined.
     */
    template <typename config1_t, typename config2_t>
    static constexpr auto is_compatible()
    {
        if constexpr (has_id_member<config1_t>(0) && has_id_member<config2_t>(0)) // needed for gcc <= 9
        {
            using config1_id_t = id_type<config1_t>;
            using config2_id_t = id_type<config2_t>;

            if constexpr (std::same_as<config1_id_t, config2_id_t>)
                return std::bool_constant<compatibility_table<config1_id_t>[as_int<config1_t>][as_int<config2_t>]>{};
            else
                return std::false_type{};
        }
        else
        {
            return std::false_type{};
        }
    }

    //!\brief Variable template that evaluates to `true` if the type has a static id member, otherwise `false`.
    template <typename config_t>
    static constexpr bool has_id = has_id_member<config_t>(0);
};
#endif // SEQAN3_WORKAROUND_GCC_PIPEABLE_CONFIG_CONCEPT

/*!\interface seqan3::detail::config_element <>
 * \brief Concept for an algorithm configuration element.
 * \ingroup algorithm
 *
 * \extends std::semiregular
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
SEQAN3_CONCEPT config_element = requires
{
    requires std::is_base_of_v<seqan3::pipeable_config_element<config_t>, config_t>;
    requires std::semiregular<config_t>;
#if SEQAN3_WORKAROUND_GCC_PIPEABLE_CONFIG_CONCEPT
    requires config_id_accessor::has_id<config_t>;
#else // ^^^ workaround / no workaround vvv
    { config_t::id };
#endif // SEQAN3_WORKAROUND_GCC_PIPEABLE_CONFIG_CONCEPT
};
//!\endcond

/*!\interface seqan3::detail::config_element_pipeable_with <>
 * \brief Concept to check if one configuration element can be combined with another configuration element.
 * \ingroup algorithm
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
SEQAN3_CONCEPT config_element_pipeable_with =
    config_element<config1_t> &&
    config_element<config2_t> &&
#if SEQAN3_WORKAROUND_GCC_PIPEABLE_CONFIG_CONCEPT
    std::same_as<config_id_accessor::id_type<config1_t>, config_id_accessor::id_type<config2_t>> &&
    decltype(config_id_accessor::is_compatible<config1_t, config2_t>())::value;
#else // ^^^ workaround / no workaround vvv
    std::same_as<std::remove_cvref_t<decltype(config1_t::id)>, std::remove_cvref_t<decltype(config2_t::id)>> &&
    compatibility_table<std::remove_cvref_t<decltype(config1_t::id)>>[static_cast<int32_t>(config1_t::id)]
                                                                     [static_cast<int32_t>(config2_t::id)];
#endif // SEQAN3_WORKAROUND_GCC_PIPEABLE_CONFIG_CONCEPT
//!\endcond

} // namespace seqan3::detail
