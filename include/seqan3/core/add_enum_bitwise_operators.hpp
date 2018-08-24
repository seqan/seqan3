// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

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
 * If this metafunction is specialised for an enum, the binary operators `&`, `|`, `^`, `~`, `&=`, `|=`, `^=` will be
 * added and behave just like for ints or unscoped enums.
 *
 * ### Example
 *
 * \snippet test/snippet/core/add_enum_bitwise_operators.cpp usage
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
