// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::add_enum_bitwise_operators.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3
{
/*!\interface seqan3::enum_bitwise_operators
 * \brief You can expect these functions on all types that overload seqan3::add_enum_bitwise_operators.
 */
/*!\name Requirements for seqan3::enum_bitwise_operators
 * \relates seqan3::enum_bitwise_operators
 * \brief You can expect these member functions.
 * \{
 * \fn operator&(t lhs, t rhs)
 * \brief Returns the binary `&` operator of lhs and rhs.
 * \param lhs First enum.
 * \param rhs Second enum.
 *
 * \returns the binary conjunction of `lhs` and `rhs`.
 *
 * \fn operator|(t lhs, t rhs)
 * \brief Returns the binary `|` operator of lhs and rhs.
 * \param lhs First enum.
 * \param rhs Second enum.
 *
 * \returns the binary disjunction of `lhs` and `rhs`.
 *
 * \fn operator^(t lhs, t rhs)
 * \brief Returns the binary `^` operator of lhs and rhs.
 * \param lhs First enum.
 * \param rhs Second enum.
 *
 * \returns the binary XOR operation on `lhs` and `rhs`.
 *
 * \fn operator~(t lhs)
 * \brief Returns the binary `~` operator of lhs.
 * \param lhs First enum.
 *
 * \returns the binary NOT operation on `lhs`.
 *
 * \fn operator&=(t & lhs, t rhs)
 * \brief Returns the binary `&=` operator of lhs and rhs.
 * \param lhs First enum.
 * \param rhs Second enum.
 *
 * \returns the binary AND assigment on `lhs`.
 *
 * \fn operator|=(t & lhs, t rhs)
 * \brief Returns the binary `|=` operator of lhs and rhs.
 * \param lhs First enum.
 * \param rhs Second enum.
 *
 * \returns the binary OR assignment on `lhs`.
 *
 * \fn operator^=(t & lhs, t rhs)
 * \brief Returns the binary `^=` operator of lhs and rhs.
 * \param lhs First enum.
 * \param rhs Second enum.
 *
 * \returns the binary XOR assignment on `lhs`.
 * \}
 */

//!\cond DEV
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
inline constexpr bool add_enum_bitwise_operators = false;

/*!\name Binary operators for scoped enums
 * \brief Perform binary operations like on ints or weak enums. These overloads are available if
 * seqan3::add_enum_bitwise_operators is defined for your type.
 * \ingroup core
 *
 * \details
 *
 * \see seqan3::add_enum_bitwise_operators
 * \{
 */
template <typename t>
constexpr t operator&(t lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    return static_cast<t>(static_cast<std::underlying_type_t<t>>(lhs) & static_cast<std::underlying_type_t<t>>(rhs));
}

template <typename t>
constexpr t operator|(t lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    return static_cast<t>(static_cast<std::underlying_type_t<t>>(lhs) | static_cast<std::underlying_type_t<t>>(rhs));
}

template <typename t>
constexpr t operator^(t lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    return static_cast<t>(static_cast<std::underlying_type_t<t>>(lhs) ^ static_cast<std::underlying_type_t<t>>(rhs));
}

template <typename t>
constexpr t operator~(t lhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    return static_cast<t>(~static_cast<std::underlying_type_t<t>>(lhs));
}

template <typename t>
constexpr t & operator&=(t & lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    lhs = lhs & rhs;
    return lhs;
}

template <typename t>
constexpr t & operator|=(t & lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    lhs = lhs | rhs;
    return lhs;
}

template <typename t>
constexpr t & operator^=(t & lhs, t rhs) noexcept
    requires std::is_enum_v<t> && add_enum_bitwise_operators<t>
{
    lhs = lhs ^ rhs;
    return lhs;
}
//!\}
//!\endcond

} // namespace seqan3
