// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the error types for maximum number of errors.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/std/concepts>

namespace seqan3::search_cfg
{

/*!\brief A strong type of underlying type `uint8_t` or `double` that represents the number or rate of total errors.
 * \ingroup search_configuration
 * \tparam value_t The underlying type.
 */
template <typename value_t>
//!\cond
    requires arithmetic<value_t>
//!\endcond
struct total : detail::strong_type<value_t, total<value_t>, detail::strong_type_skill::convert>
{
    using detail::strong_type<value_t, total<value_t>, detail::strong_type_skill::convert>::strong_type;

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr std::integral_constant<uint8_t, 0> _id{};
};

/*!\name Template argument type deduction guides
 * \relates seqan3::search_cfg::total
 * \{
 */
//! \brief Deduces to `uint8_t` for all types modelling std::integral.
template <std::integral value_t>
total(value_t) -> total<uint8_t>;

//! \brief Deduces to `double` for all types modelling seqan3::floating_point.
template <typename value_t>
//!\cond
    requires floating_point<value_t>
//!\endcond
total(value_t) -> total<double>;
//!\}

/*!\brief A strong type of underlying type `uint8_t` or `double` that represents the number or rate of
 *        substitutions.
 * \ingroup search_configuration
 * \tparam value_t The underlying type.
 */
template <typename value_t>
//!\cond
    requires arithmetic<value_t>
//!\endcond
struct substitution : detail::strong_type<value_t, substitution<value_t>, detail::strong_type_skill::convert>
{
    using detail::strong_type<value_t, substitution<value_t>, detail::strong_type_skill::convert>::strong_type;

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr std::integral_constant<uint8_t, 1> _id{};
};

/*!\name Template argument type deduction guides
 * \relates seqan3::search_cfg::substitution
 * \{
 */
//! \brief Deduces to `uint8_t` for all types modelling std::integral.
template <std::integral value_t>
substitution(value_t) -> substitution<uint8_t>;

//! \brief Deduces to `double` for all types modelling seqan3::floating_point.
template <typename value_t>
//!\cond
    requires floating_point<value_t>
//!\endcond
substitution(value_t) -> substitution<double>;
//!\}

/*!\brief A strong type of underlying type `uint8_t` or `double` that represents the number or rate of insertions.
 * \ingroup search_configuration
 * \tparam value_t The underlying type.
 */
template <typename value_t>
//!\cond
    requires arithmetic<value_t>
//!\endcond
struct insertion : detail::strong_type<value_t, insertion<value_t>, detail::strong_type_skill::convert>
{
    using detail::strong_type<value_t, insertion<value_t>, detail::strong_type_skill::convert>::strong_type;

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr std::integral_constant<uint8_t, 2> _id{};
};

/*!\name Template argument type deduction guides
 * \relates seqan3::search_cfg::total
 * \{
 */
//! \brief Deduces to `uint8_t` for all types modelling std::integral.
template <std::integral value_t>
insertion(value_t) -> insertion<uint8_t>;

//! \brief Deduces to `double` for all types modelling seqan3::floating_point.
template <typename value_t>
//!\cond
    requires floating_point<value_t>
//!\endcond
insertion(value_t) -> insertion<double>;
//!\}

/*!\brief A strong type of underlying type `uint8_t` or `double` that represents the number or rate of deletions.
 * \ingroup search_configuration
 * \tparam value_t The underlying type.
 */
template <typename value_t>
//!\cond
    requires arithmetic<value_t>
//!\endcond
struct deletion : detail::strong_type<value_t, deletion<value_t>, detail::strong_type_skill::convert>
{
    using detail::strong_type<value_t, deletion<value_t>, detail::strong_type_skill::convert>::strong_type;

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr std::integral_constant<uint8_t, 3> _id{};
};

/*!\name Template argument type deduction guides
 * \relates seqan3::search_cfg::deletion
 * \{
 */
//! \brief Deduces to `uint8_t` for integral types.
template <std::integral value_t>
deletion(value_t) -> deletion<uint8_t>;

//! \brief Deduces to `double` for all types modelling seqan3::floating_point.
template <typename value_t>
//!\cond
    requires floating_point<value_t>
//!\endcond
deletion(value_t) -> deletion<double>;
//!\}

} // namespace seqan3::search_cfg
