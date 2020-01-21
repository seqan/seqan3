// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::add_enum_bitwise_operators.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief Set to true for a scoped enum to have binary operators overloaded.
 * \ingroup core
 *
 * \details
 *
 * If this type trait is specialised for an enum, the binary operators `&`, `|`, `^`, `~`, `&=`, `|=`, `^=` will be
 * added and behave just like for ints or unscoped enums.
 *
 * ### Example
 *
 * \include test/snippet/core/add_enum_bitwise_operators.cpp
 */
template <typename t>
constexpr bool add_enum_bitwise_operators = false;

/*!\name Binary operators for scoped enums
 * \brief Perform binary operations like on ints or weak enums. These overloads are available if
 * seqan3::add_enum_bitwise_operators is defined for your type.
 * \ingroup core
 * \see seqan3::add_enum_bitwise_operators
 * \{
 */
template <typename t>
constexpr t operator& (t lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    return static_cast<t>(static_cast<std::underlying_type_t<t>>(lhs) & static_cast<std::underlying_type_t<t>>(rhs));
}

template <typename t>
constexpr t operator| (t lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    return static_cast<t>(static_cast<std::underlying_type_t<t>>(lhs) | static_cast<std::underlying_type_t<t>>(rhs));
}

template <typename t>
constexpr t operator^ (t lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    return static_cast<t>(static_cast<std::underlying_type_t<t>>(lhs) ^ static_cast<std::underlying_type_t<t>>(rhs));
}

template <typename t>
constexpr t operator~ (t lhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    return static_cast<t>(~static_cast<std::underlying_type_t<t>>(lhs));
}

template <typename t>
constexpr t & operator&= (t & lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    lhs = lhs & rhs;
    return lhs;
}

template <typename t>
constexpr t & operator|= (t & lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    lhs = lhs | rhs;
    return lhs;
}

template <typename t>
constexpr t & operator^= (t & lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    lhs = lhs ^ rhs;
    return lhs;
}
//!\}

} // namespace seqan3
