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

#include <functional>
#include <type_traits>

#include <meta/meta.hpp>

#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/std/concepts>

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
inline constexpr std::array<std::array<void, 0>, 0> compatibility_table;

// ----------------------------------------------------------------------------
// Concept config_element_specialisation
// ----------------------------------------------------------------------------

#ifdef SEQAN3_WORKAROUND_GCC_PIPEABLE_CONFIG_CONCEPT
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
    //!\brief Variable template that evaluates to `true` if the type has a static id member, otherwise `false`.
    template <typename config_t>
    static constexpr bool has_id = has_id_member<config_t>(0);
};
#endif // SEQAN3_WORKAROUND_GCC_PIPEABLE_CONFIG_CONCEPT

/*!\interface seqan3::detail::config_element_specialisation <>
 * \brief Concept for an algorithm configuration element.
 * \ingroup algorithm
 *
 * \extends std::semiregular
 * \implements seqan3::pipeable_config_element
 */

/*!\name Requirements for seqan3::detail::config_element_specialisation
 * \relates seqan3::detail::config_element_specialisation
 * \brief   You can expect this member on all types that satisfy seqan3::detail::config_element_specialisation.
 * \{
 */
/*!\var id
 * \brief Algorithm specific static id used for internal validation checks.
 */
//!\}
//!\cond
template <typename config_t>
SEQAN3_CONCEPT config_element_specialisation = requires
{
    requires std::is_base_of_v<seqan3::pipeable_config_element<config_t>, config_t>;
    requires std::semiregular<config_t>;
#ifdef SEQAN3_WORKAROUND_GCC_PIPEABLE_CONFIG_CONCEPT
    requires config_id_accessor::has_id<config_t>;
#else // ^^^ workaround / no workaround vvv
    { config_t::id };
#endif // SEQAN3_WORKAROUND_GCC_PIPEABLE_CONFIG_CONCEPT
};
//!\endcond

} // namespace seqan3::detail
