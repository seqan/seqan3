// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concepts for core language types and relations that don't have concepts in C++20 (yet).
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/platform.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\addtogroup concept
 * \{
 */

/*!\interface   seqan3::detail::WeaklyEqualityComparableByMembersWith <>
 * \brief       Like std::detail::WeaklyEqualityComparableWith, but considers only member operators of the LHS.
 */
//!\cond
template <typename lhs_t, typename rhs_t>
SEQAN3_CONCEPT WeaklyEqualityComparableByMembersWith = requires (lhs_t const & lhs, rhs_t const & rhs)
{
    lhs.operator==(rhs); std::Boolean<decltype(lhs.operator==(rhs))>;
    lhs.operator!=(rhs); std::Boolean<decltype(lhs.operator!=(rhs))>;
};
//!\endcond
/*!\interface   seqan3::detail::WeaklyOrderedByMembersWith <>
 * \brief       Like seqan3::WeaklyOrderedWith, but considers only member operators of the LHS.
 */
//!\cond
template <typename lhs_t, typename rhs_t>
SEQAN3_CONCEPT WeaklyOrderedByMembersWith = requires (lhs_t const & lhs, rhs_t const & rhs)
{
    lhs.operator< (rhs); std::Boolean<decltype(lhs.operator< (rhs))>;
    lhs.operator> (rhs); std::Boolean<decltype(lhs.operator> (rhs))>;
    lhs.operator<=(rhs); std::Boolean<decltype(lhs.operator<=(rhs))>;
    lhs.operator>=(rhs); std::Boolean<decltype(lhs.operator>=(rhs))>;
};
//!\endcond

/*!\interface   seqan3::detail::ConvertibleToByMember <>
 * \brief       Like seqan3::ImplicitlyConvertibleTo, but only considers member operators of the source type.
 */
//!\cond
template <typename source_t, typename target_t>
SEQAN3_CONCEPT ConvertibleToByMember = requires (source_t s)
{
    { s.operator target_t() } -> target_t;
};
//!\endcond
//!\}

} // seqan3::detail

namespace seqan3
{

/*!\addtogroup concept
 * \{
 */

/*!\interface   seqan3::WeaklyOrderedWith <>
 * \tparam t1   The first type to compare.
 * \tparam t2   The second type to compare.
 * \brief       Requires the two operands to be comparable with `==` and `!=` in both directions.
 * \sa          http://en.cppreference.com/w/cpp/experimental/ranges/concepts/WeaklyEqualityComparableWith
 */
//!\cond
template <typename t1, typename t2>
SEQAN3_CONCEPT WeaklyOrderedWith = requires (std::remove_reference_t<t1> const & v1,
                                             std::remove_reference_t<t2> const & v2)
{
    { v1 <  v2 } -> bool &&;
    { v1 <= v2 } -> bool &&;
    { v2 >  v1 } -> bool &&;
    { v2 >= v1 } -> bool &&;
};
//!\endcond

/*!\interface   seqan3::ImplicitlyConvertibleTo <>
 * \brief       Resolves to `std::ranges::ImplicitlyConvertibleTo<type1, type2>()`.
 */
//!\cond
template <typename t, typename u>
SEQAN3_CONCEPT ImplicitlyConvertibleTo = std::is_convertible_v<t, u>;
//!\endcond

/*!\interface   seqan3::ExplicitlyConvertibleTo <>
 * \brief       Resolves to `std::ranges::ExplicitlyConvertibleTo<type1, type2>()`.
 */
//!\cond
template <typename t, typename u>
SEQAN3_CONCEPT ExplicitlyConvertibleTo = requires (t vt) { { static_cast<u>(vt)}; };
//!\endcond

/*!\interface   seqan3::Arithmetic <>
 * \brief       A type that satisfies std::is_arithmetic_v<t>.
 * \sa          http://en.cppreference.com/w/cpp/types/is_arithmetic
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT Arithmetic = std::is_arithmetic_v<t>;
//!\endcond

/*!\interface   seqan3::FloatingPoint <>
 * \extends     seqan3::Arithmetic
 * \brief       An arithmetic type that also satisfies std::is_floating_point_v<t>.
 * \sa          http://en.cppreference.com/w/cpp/types/is_floating_point
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT FloatingPoint = Arithmetic<t> && std::is_floating_point_v<t>;
//!\endcond

/*!\interface   seqan3::Char <>
 * \extends     std::Integral
 * \brief       This concept encompasses exactly the types `char`, `signed char`, `unsigned char`, `wchar_t`,
 *              `char16_t` and `char32_t`.
 */
//!\cond

template <typename t>
SEQAN3_CONCEPT Char = std::Integral<t> &&
                      (std::Same<t, char> || std::Same<t, unsigned char> || std::Same<t, signed char> ||
#ifdef __cpp_char8_t
                       std::Same<t, char8_t> ||
#endif
                       std::Same<t, char16_t> || std::Same<t, char32_t> || std::Same<t, wchar_t>);
//!\endcond

/*!\interface   seqan3::TriviallyDestructible <>
 * \extends     std::Destructible
 * \brief       A type that satisfies std::is_trivially_destructible_v<t>.
 * \sa          http://en.cppreference.com/w/cpp/types/is_destructible
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT TriviallyDestructible = std::Destructible<t> && std::is_trivially_destructible_v<t>;
//!\endcond

/*!\interface   seqan3::TriviallyCopyable
 * \brief       A type that satisfies std::is_trivially_copyable_v<t>.
 * \extends     std::Copyable
 * \sa          http://en.cppreference.com/w/cpp/types/is_trivially_copyable
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT TriviallyCopyable = std::Copyable<t> && std::is_trivially_copyable_v<t>;
//!\endcond

/*!\interface   seqan3::Trivial
 * \brief       A type that satisfies seqan3::TriviallyCopyable and seqan3::TriviallyDestructible.
 * \extends     seqan3::TriviallyCopyable
 * \extends     seqan3::TriviallyDestructible
 * \sa          http://en.cppreference.com/w/cpp/experimental/ranges/concepts/Copyable
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT Trivial = TriviallyCopyable<t> && TriviallyDestructible<t> && std::is_trivial_v<t>;
//!\endcond

/*!\interface   seqan3::StandardLayout
 * \brief       Resolves to std::is_standard_layout_v<t>.
 * \sa          https://en.cppreference.com/w/cpp/named_req/StandardLayoutType
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT StandardLayout = std::is_standard_layout_v<t>;
//!\endcond

/*!\interface   seqan3::WeaklyAssignable
 * \brief       Resolves to std::is_assignable_v<t>.
 * \sa          https://en.cppreference.com/w/cpp/types/is_assignable
 *
 * \details
 *
 * Note that this requires less than std::Assignable, it simply tests if the expression
 * `std::declval<T>() = std::declval<U>()` is well-formed.
 */
//!\cond
template <typename t, typename u>
SEQAN3_CONCEPT WeaklyAssignable = std::is_assignable_v<t, u>;
//!\endcond

}  // namespace seqan3
