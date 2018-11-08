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
 * \brief Provides utilities for modifying characters.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <cmath>

#include <seqan3/core/detail/int_types.hpp>
#include <seqan3/core/concept/core_language.hpp>

namespace seqan3::detail
{

//!\brief Auxiliary table for seqan3::to_lower.
template <typename char_type>
inline std::array<char_type, detail::size_in_values_v<char_type>> constexpr to_lower_table
{
    [] () constexpr
    {
        std::array<char_type, detail::size_in_values_v<char_type>> ret{};

        for (size_t i = 0; i < detail::size_in_values_v<char_type>; ++i)
            ret[i] = i;

        for (size_t i = char_type{'A'}; i <= char_type{'Z'}; ++i)
            ret[i] = ret[i] - char_type{'A'} + char_type{'a'};

        return ret;
    } ()
};

//!\brief Auxiliary table for seqan3::to_upper.
template <typename char_type>
inline std::array<char_type, detail::size_in_values_v<char_type>> constexpr to_upper_table
{
    [] () constexpr
    {
        std::array<char_type, detail::size_in_values_v<char_type>> ret{};

        for (size_t i = 0; i < detail::size_in_values_v<char_type>; ++i)
            ret[i] = i;

        for (size_t i = char_type{'a'}; i <= char_type{'z'}; ++i)
            ret[i] = ret[i] - char_type{'a'} + char_type{'A'};

        return ret;
    } ()
};

} // namespace seqan3::detail

namespace seqan3
{

/*!\name Operations on characters
 * \ingroup stream
 * \{
 */

/*!\brief Converts 'A'-'Z' to 'a'-'z' respectively; other characters are returned as is.
 * \tparam char_type Type of the parameter; must model seqan3::char_concept.
 * \param c The parameter.
 * \returns The character converted to lower case.
 *
 * \details
 *
 * In contrast to std::tolower this function is independent of locale and can be evaluated in a `constexpr` context.
 */
template <char_concept char_type>
constexpr char_type to_lower(char_type const c) noexcept
{
    using u_t = std::make_unsigned_t<char_type>;
    return detail::to_lower_table<char_type>[static_cast<u_t>(c)];
}

/*!\brief Converts 'a'-'z' to 'A'-'Z' respectively; other characters are returned as is.
 * \tparam char_type Type of the parameter; must model seqan3::char_concept.
 * \param c The parameter.
 * \returns The character converted to upper case.
 *
 * \details
 *
 * In contrast to std::to_upper this function is independent of locale and can be evaluated in a `constexpr` context.
 */
template <char_concept char_type>
constexpr char_type to_upper(char_type const c) noexcept
{
    using u_t = std::make_unsigned_t<char_type>;
    return detail::to_upper_table<char_type>[static_cast<u_t>(c)];
}
//!\}

} // namespace seqan3
