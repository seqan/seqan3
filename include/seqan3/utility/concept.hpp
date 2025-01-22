// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides concepts that do not have equivalents in C++20.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <seqan3/core/platform.hpp>

/*!\defgroup utility_concept Concept
 * \brief Provides various general purpose concepts.
 * \ingroup utility
 * \see utility
 */

namespace seqan3
{
/*!\interface seqan3::implicitly_convertible_to <>
 * \ingroup utility_concept
 * \brief Checks whether `from` can be implicityly converted to `to`.
 *
 * \details
 *
 * The STL concept `convertible_to` checks for implicit and explicit convertability.
 *
 * \noapi
 */
//!\cond
template <typename from, typename to>
concept implicitly_convertible_to = std::is_convertible_v<from, to>;
//!\endcond

/*!\interface seqan3::explicitly_convertible_to <>
 * \ingroup utility_concept
 * \brief Checks whether `from` can be explicitly converted to `to`.
 *
 * \details
 *
 * The STL concept `convertible_to` checks for implicit and explicit convertability.
 *
 * \noapi
 */
//!\cond
template <typename from, typename to>
concept explicitly_convertible_to = requires { static_cast<to>(std::declval<from>()); };
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
concept builtin_character = std::integral<t>
                         && (std::same_as<t, char> || std::same_as<t, unsigned char> || std::same_as<t, signed char> ||
#ifdef __cpp_char8_t
                             std::same_as<t, char8_t> ||
#endif
                             std::same_as<t, char16_t> || std::same_as<t, char32_t> || std::same_as<t, wchar_t>);
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
concept trivial = trivially_copyable<t> && trivially_destructible<t> && std::is_trivially_default_constructible_v<t>;
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
} // namespace seqan3
