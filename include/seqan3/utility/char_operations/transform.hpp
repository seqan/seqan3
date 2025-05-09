// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides utilities for modifying characters.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <array>

#include <seqan3/utility/concept.hpp>
#include <seqan3/utility/detail/integer_traits.hpp>

namespace seqan3::detail
{

//!\brief Auxiliary table for seqan3::to_lower.
//!\ingroup utility_char_operations
template <typename char_type>
inline constexpr std::array<char_type, detail::size_in_values_v<char_type>> to_lower_table{
    []() constexpr
    {
        std::array<char_type, detail::size_in_values_v<char_type>> ret{};

        for (size_t i = 0; i < detail::size_in_values_v<char_type>; ++i)
            ret[i] = i;

        for (size_t i = char_type{'A'}; i <= char_type{'Z'}; ++i)
            ret[i] = ret[i] - char_type{'A'} + char_type{'a'};

        return ret;
    }()};

//!\brief Auxiliary table for seqan3::to_upper.
//!\ingroup utility_char_operations
template <typename char_type>
inline constexpr std::array<char_type, detail::size_in_values_v<char_type>> to_upper_table{
    []() constexpr
    {
        std::array<char_type, detail::size_in_values_v<char_type>> ret{};

        for (size_t i = 0; i < detail::size_in_values_v<char_type>; ++i)
            ret[i] = i;

        for (size_t i = char_type{'a'}; i <= char_type{'z'}; ++i)
            ret[i] = ret[i] - char_type{'a'} + char_type{'A'};

        return ret;
    }()};

} // namespace seqan3::detail

namespace seqan3
{

/*!\name Operations on characters
 * \ingroup utility_char_operations
 * \{
 */

/*!\brief Converts 'A'-'Z' to 'a'-'z' respectively; other characters are returned as is.
 * \tparam char_type Type of the parameter; must model seqan3::builtin_character.
 * \param c The parameter.
 * \returns The character converted to lower case.
 *
 * \details
 *
 * In contrast to std::tolower this function is independent of locale and can be evaluated in a `constexpr` context.
 */
template <builtin_character char_type>
constexpr char_type to_lower(char_type const c) noexcept
{
    using u_t = std::make_unsigned_t<char_type>;
    return detail::to_lower_table<char_type>[static_cast<u_t>(c)];
}

/*!\brief Converts 'a'-'z' to 'A'-'Z' respectively; other characters are returned as is.
 * \tparam char_type Type of the parameter; must model seqan3::builtin_character.
 * \param c The parameter.
 * \returns The character converted to upper case.
 *
 * \details
 *
 * In contrast to std::to_upper this function is independent of locale and can be evaluated in a `constexpr` context.
 */
template <builtin_character char_type>
constexpr char_type to_upper(char_type const c) noexcept
{
    using u_t = std::make_unsigned_t<char_type>;
    return detail::to_upper_table<char_type>[static_cast<u_t>(c)];
}
//!\}

} // namespace seqan3
