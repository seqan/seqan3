// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concepts for core language types and relations that don't have concepts in C++20 (yet).
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concepts>
#include <type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{
/*!\interface   seqan3::detail::weakly_equality_comparable_with <>
 * \ingroup utility_concept
 * \tparam t1   The first type to compare.
 * \tparam t2   The second type to compare.
 * \brief       Requires the two operands to be comparable with `==` and `!=` in both directions.
 */
//!\cond
template <class T, class U>
concept weakly_equality_comparable_with =
    requires(std::remove_reference_t<T> const & t,
             std::remove_reference_t<U> const & u)
    {
        requires std::convertible_to<decltype(t == u), bool>;
        requires std::convertible_to<decltype(t != u), bool>;
        requires std::convertible_to<decltype(u == t), bool>;
        requires std::convertible_to<decltype(u != t), bool>;
    };
//!\endcond

//!\brief Binary type trait that behaves like the seqan3::detail::weakly_equality_comparable_with concept.
//!\ingroup utility_concept
template <typename lhs_t, typename rhs_t>
struct weakly_equality_comparable_with_trait :
    std::integral_constant<bool, weakly_equality_comparable_with<lhs_t, rhs_t>>
{};

/*!\interface   seqan3::detail::weakly_ordered_with <>
 * \ingroup utility_concept
 * \tparam t1   The first type to compare.
 * \tparam t2   The second type to compare.
 * \brief       Requires the two operands to be comparable with `<`, `<=`, `>` and `>=` in both directions.
 */
//!\cond
template <typename t1, typename t2>
concept weakly_ordered_with = requires (std::remove_reference_t<t1> const & v1,
                                        std::remove_reference_t<t2> const & v2)
{
    requires std::convertible_to<decltype(v1 <  v2), bool>;
    requires std::convertible_to<decltype(v1 <= v2), bool>;
    requires std::convertible_to<decltype(v1 >  v2), bool>;
    requires std::convertible_to<decltype(v1 >= v2), bool>;

    requires std::convertible_to<decltype(v2 <  v1), bool>;
    requires std::convertible_to<decltype(v2 <= v1), bool>;
    requires std::convertible_to<decltype(v2 >  v1), bool>;
    requires std::convertible_to<decltype(v2 >= v1), bool>;
};
//!\endcond

//!\brief Binary type trait that behaves like the seqan3::detail::weakly_ordered_with concept.
//!\ingroup utility_concept
template <typename lhs_t, typename rhs_t>
struct weakly_ordered_with_trait : std::integral_constant<bool, weakly_ordered_with<lhs_t, rhs_t>>
{};

} // seqan3::detail

namespace seqan3
{
/*!\interface seqan3::implicitly_convertible_to <>
 * \ingroup utility_concept
 * \brief Resolves to `std::ranges::implicitly_convertible_to<type1, type2>()`.
 *
 * \details
 *
 * \noapi
 */
//!\cond
template <typename t, typename u>
concept implicitly_convertible_to = std::is_convertible_v<t, u>;
//!\endcond

/*!\interface seqan3::explicitly_convertible_to <>
 * \ingroup utility_concept
 * \brief Resolves to `std::ranges::explicitly_convertible_to<type1, type2>()`.
 *
 * \details
 *
 * \noapi
 */
//!\cond
template <typename t, typename u>
concept explicitly_convertible_to = requires (t vt) { { static_cast<u>(vt)}; };
//!\endcond

/*!\interface seqan3::arithmetic <>
 * \ingroup utility_concept
 * \brief A type that satisfies std::is_arithmetic_v<t>.
 *
 * \details
 *
 * \sa https://en.cppreference.com/w/cpp/types/is_arithmetic
 *
 * \noapi
 */
//!\cond
template <typename t>
concept arithmetic = std::is_arithmetic_v<t>;
//!\endcond

/*!\interface seqan3::floating_point <>
 * \ingroup utility_concept
 * \extends seqan3::arithmetic
 * \brief An arithmetic type that also satisfies std::is_floating_point_v<t>.
 *
 * \details
 *
 * \sa https://en.cppreference.com/w/cpp/types/is_floating_point
 *
 * \noapi
 */
//!\cond
template <typename t>
concept floating_point = arithmetic<t> && std::is_floating_point_v<t>;
//!\endcond

/*!\interface seqan3::builtin_character <>
 * \ingroup utility_concept
 * \extends std::integral
 * \brief This concept encompasses exactly the types `char`, `signed char`, `unsigned char`, `wchar_t`,
 *        `char16_t` and `char32_t`.
 *
 * \details
 *
 * \noapi
 */
//!\cond

template <typename t>
concept builtin_character = std::integral<t> &&
                            (std::same_as<t, char> ||
                             std::same_as<t, unsigned char> ||
                             std::same_as<t, signed char> ||
#ifdef __cpp_char8_t
                             std::same_as<t, char8_t> ||
#endif
                             std::same_as<t, char16_t> ||
                             std::same_as<t, char32_t> ||
                             std::same_as<t, wchar_t>);
//!\endcond

/*!\interface seqan3::trivially_destructible <>
 * \ingroup utility_concept
 * \extends std::destructible
 * \brief A type that satisfies std::is_trivially_destructible_v<t>.
 *
 * \details
 *
 * \sa https://en.cppreference.com/w/cpp/types/is_destructible
 *
 * \noapi
 */
//!\cond
template <typename t>
concept trivially_destructible = std::destructible<t> && std::is_trivially_destructible_v<t>;
//!\endcond

/*!\interface seqan3::trivially_copyable
 * \ingroup utility_concept
 * \brief A type that satisfies std::is_trivially_copyable_v<t>.
 * \extends std::copyable
 *
 * \details
 *
 * \sa https://en.cppreference.com/w/cpp/types/is_trivially_copyable
 *
 * \noapi
 */
//!\cond
template <typename t>
concept trivially_copyable = std::copyable<t> && std::is_trivially_copyable_v<t>;
//!\endcond

/*!\interface seqan3::trivial
 * \ingroup utility_concept
 * \brief A type that satisfies seqan3::trivially_copyable and seqan3::trivially_destructible.
 * \extends seqan3::trivially_copyable
 * \extends seqan3::trivially_destructible
 *
 * \details
 *
 * \sa https://en.cppreference.com/w/cpp/types/is_trivial
 *
 * \noapi
 */
//!\cond
template <typename t>
concept trivial = trivially_copyable<t> && trivially_destructible<t> && std::is_trivial_v<t>;
//!\endcond

/*!\interface seqan3::standard_layout
 * \ingroup utility_concept
 * \brief Resolves to std::is_standard_layout_v<t>.
 *
 * \details
 *
 * \sa https://en.cppreference.com/w/cpp/types/is_standard_layout
 *
 * \noapi
 */
//!\cond
template <typename t>
concept standard_layout = std::is_standard_layout_v<t>;
//!\endcond

/*!\interface seqan3::weakly_assignable_from
 * \ingroup utility_concept
 * \brief Resolves to std::is_assignable_v<t>.
 *
 * \details
 *
 * \note This requires less than std::assignable_from, it simply tests if the expression `std::declval<T>() =
 *       std::declval<U>()` is well-formed.
 *
 * \sa https://en.cppreference.com/w/cpp/types/is_assignable
 *
 * \noapi
 */
//!\cond
template <typename t, typename u>
concept weakly_assignable_from = std::is_assignable_v<t, u>;
//!\endcond
}  // namespace seqan3
