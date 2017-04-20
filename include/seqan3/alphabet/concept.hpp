// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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

/*!\file alphabet/concept.hpp
 * \ingroup alphabet
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Core alphabet concept.
 */

#pragma once

#include <iostream>
#include <string>

namespace seqan3
{

// ------------------------------------------------------------------
// type metafunctions operator
// ------------------------------------------------------------------

//! Type metafunction that returns the `char_type` defined inside an alphabet type.
template <typename alphabet_type>
    requires requires (alphabet_type c) { typename alphabet_type::char_type; }
struct underlying_char
{
    using type = typename alphabet_type::char_type;
};

//! Shortcut for @link underlying_char @endlink
template <typename alphabet_type>
using underlying_char_t = typename underlying_char<alphabet_type>::type;

//! Type metafunction that returns the `rank_type` defined inside an alphabet type.
template <typename alphabet_type>
    requires requires (alphabet_type c) { typename alphabet_type::rank_type; }
struct underlying_rank
{
    using type = typename alphabet_type::rank_type;
};

//! Shortcut for @link underlying_rank @endlink
template <typename alphabet_type>
using underlying_rank_t = typename underlying_rank<alphabet_type>::type;

// ------------------------------------------------------------------
// value metafunctions
// ------------------------------------------------------------------

//! Value metafunction that returns the `value_size` defined inside an alphabet type.
template <typename alphabet_type>
    requires requires (alphabet_type c) { alphabet_type::value_size; }
struct alphabet_size
{
    static constexpr underlying_rank_t<alphabet_type> value = alphabet_type::value_size;
};

//! Shortcut for @link alphabet_size @endlink
template <typename alphabet_type>
constexpr underlying_rank_t<alphabet_type> alphabet_size_v = alphabet_size<alphabet_type>::value;

// ------------------------------------------------------------------
// free functions
// ------------------------------------------------------------------

//!\publicsection
//!@name Wrapper functions to make alphabet members "globally" visible
//!@{

//! Free function that calls `.to_char()` on the argument
template <typename alphabet_type>
constexpr underlying_char_t<alphabet_type> to_char(alphabet_type const & c)
    requires requires (alphabet_type c) { { c.to_char() } -> underlying_char_t<alphabet_type>; }
{
    return c.to_char();
}

//! Free function that calls `.to_rank()` on the argument
template <typename alphabet_type>
constexpr underlying_rank_t<alphabet_type> to_rank(alphabet_type const & c)
    requires requires (alphabet_type c) { { c.to_rank() } -> underlying_rank_t<alphabet_type>; }
{
    return c.to_rank();
}

//! Free function that calls `.assign_char(in)` on the first argument
template <typename alphabet_type>
constexpr alphabet_type & assign_char(alphabet_type & c, char const in)
    requires requires (alphabet_type c) { { c.assign_char(char{0}) } -> alphabet_type &; }
{
    return c.assign_char(in);
}

//! Free function that calls `.assign_rank(in)` on the first argument
template <typename alphabet_type>
constexpr alphabet_type & assign_rank(alphabet_type & c, underlying_rank_t<alphabet_type> const in)
    requires requires (alphabet_type c) { { c.assign_rank(uint8_t{0}) } -> alphabet_type &; }
{
    return c.assign_rank(in);
}

//! Free ostream operator that delegates to `c.to_char()`
template <typename alphabet_type>
std::ostream& operator<<(std::ostream & os, alphabet_type const & c)
    requires requires (alphabet_type c) { { c.to_char() } -> underlying_char_t<alphabet_type>; }
{
    os << c.to_char();
    return os;
}

//!@}
// ------------------------------------------------------------------
// alphabet concept
// ------------------------------------------------------------------

/*!\var concept bool alphabet_concept
 * \brief A concept for container and string alphabets
 * \privatesection
 */

template <typename t>
concept bool alphabet_concept = requires (t t1, t t2)
{
    // StL concepts
    requires std::is_pod_v<t> == true;
    requires std::is_swappable_v<t> == true;

    // static data members
    alphabet_size<t>::value;

    // conversion from/to char
    { to_char(t1)     } -> underlying_char_t<t>;
    { to_rank(t1) } -> underlying_rank_t<t>;

    { assign_char(t1, 0)     } -> t;
    { assign_rank(t1, 0) } -> t;

    { std::cout << t1 };

    // required comparison operators
    { t1 == t2 } -> bool;
    { t1 != t2 } -> bool;
    { t1 <  t2 } -> bool;
    { t1 >  t2 } -> bool;
    { t1 <= t2 } -> bool;
    { t1 >= t2 } -> bool;
};

//TODO serialization

} // namespace seqan3
