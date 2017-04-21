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
 * \brief Core alphabet concept and free function/metafunction wrappers.
 */

#pragma once

#include <iostream>
#include <string>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

// ============================================================================
// Free/Global interface wrappers
// ============================================================================

//!\addtogroup alphabet
//!\{

// ------------------------------------------------------------------
// type metafunctions
// ------------------------------------------------------------------

//!\brief Type metafunction that returns the `char_type` defined inside an alphabet type.
//!\tparam alphabet_type Must provide a `char_type` member type.
template <typename alphabet_type>
//!\cond
    requires requires (alphabet_type c) { typename alphabet_type::char_type; }
//!\endcond
struct underlying_char
{
    //!\brief The alphabet's char_type.
    using type = typename alphabet_type::char_type;
};

//!\brief Shortcut for seqan3::underlying_char
template <typename alphabet_type>
//!\relates seqan3::underlying_char
using underlying_char_t = typename underlying_char<alphabet_type>::type;

//!\brief Type metafunction that returns the `rank_type` defined inside an alphabet type.
//!\tparam alphabet_type Must provide a `rank_type` member type.
template <typename alphabet_type>
//!\cond
    requires requires (alphabet_type c) { typename alphabet_type::rank_type; }
//!\endcond
struct underlying_rank
{
    //!\brief The alphabet's rank_type.
    using type = typename alphabet_type::rank_type;
};

//!\brief Shortcut for seqan3::underlying_rank
//!\relates seqan3::underlying_rank
template <typename alphabet_type>
using underlying_rank_t = typename underlying_rank<alphabet_type>::type;

// ------------------------------------------------------------------
// value metafunctions
// ------------------------------------------------------------------

//!\brief Value metafunction that returns the `value_size` defined inside an alphabet type.
//!\tparam alphabet_type Must provide a `value_size` static member.
template <typename alphabet_type>
//!\cond
    requires requires (alphabet_type c) { alphabet_type::value_size; }
//!\endcond
struct alphabet_size
{
    //!\brief The alphabet's size.
    static constexpr underlying_rank_t<alphabet_type> value = alphabet_type::value_size;
};

//!\brief Shortcut for seqan3::alphabet_size
//!\relates seqan3::alphabet_size
template <typename alphabet_type>
constexpr underlying_rank_t<alphabet_type> alphabet_size_v = alphabet_size<alphabet_type>::value;

// ------------------------------------------------------------------
// free functions
// ------------------------------------------------------------------

/*!\name Free function wrappers for alphabet member functions
 * \brief For alphabets that implement needed functions as members, make them "globally" available.
 * \{
 */
/*!\brief Free function wrapper that calls `.to_char()` on the argument.
 * \tparam alphabet_type Must provide a `.to_char()` member function.
 * \param c The alphabet letter that you wish to convert to char.
 * \returns The letter's value in the alphabet's rank type (usually `char`).
 */
template <typename alphabet_type>
constexpr underlying_char_t<alphabet_type> to_char(alphabet_type const c)
    requires requires (alphabet_type c) { { c.to_char() } -> underlying_char_t<alphabet_type>; }
{
    return c.to_char();
}

/*!\brief Ostream operator that delegates to `c.to_char()`.
 * \tparam alphabet_type Must provide a `.to_char()` member function.
 * \param os The output stream you are printing to.
 * \param c The alphabet letter that you wish to convert to char.
 * \returns A reference to the output stream.
 */
template <typename alphabet_type>
std::ostream & operator<<(std::ostream & os, alphabet_type const c)
//!\cond
    requires requires (alphabet_type c) { { c.to_char() } -> underlying_char_t<alphabet_type>; }
//!\endcond
{
    os << c.to_char();
    return os;
}

/*!\brief Free function wrapper that calls `.to_rank()` on the argument.
 * \tparam alphabet_type Must provide a `.to_rank()` member function.
 * \param c The alphabet letter that you wish to convert to rank.
 * \returns The letter's value in the alphabet's rank type (usually a `uint*_t`).
 */
template <typename alphabet_type>
constexpr underlying_rank_t<alphabet_type> to_rank(alphabet_type const c)
    requires requires (alphabet_type c) { { c.to_rank() } -> underlying_rank_t<alphabet_type>; }
{
    return c.to_rank();
}

/*!\brief Free function wrapper that calls `.assign_char(in)` on the first argument.
 * \tparam alphabet_type Must provide an `.assign_char()` member function.
 * \param c The alphabet letter that you wish to assign to.
 * \param in The `char` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename alphabet_type>
constexpr alphabet_type & assign_char(alphabet_type & c, underlying_char_t<alphabet_type> const in)
    requires requires (alphabet_type c) { { c.assign_char(char{0}) } -> alphabet_type &; }
{
    return c.assign_char(in);
}

/*!\brief Free function wrapper that calls `.assign_char(in)` on the first argument.
 * \tparam alphabet_type Must provide an `.assign_char()` member function.
 * \param c An alphabet letter temporary.
 * \param in The `char` value you wish to assign.
 * \returns The assignment result as a temporary.
 * \details
 * Use this e.g. to newly create alphabet letters from char:
 * ~~~{.cpp}
 * auto l = assign_char(dna5{}, 'G');  // l is of type dna5
 * ~~~
 */
template <typename alphabet_type>
constexpr alphabet_type && assign_char(alphabet_type && c, underlying_char_t<alphabet_type> const in)
    requires requires (alphabet_type c) { { c.assign_char(char{0}) } -> alphabet_type &; }
{
    return std::move(c.assign_char(in));
}

/*!\brief Free function wrapper that calls `.assign_rank(in)` on the first argument.
 * \tparam alphabet_type Must provide an `.assign_rank()` member function.
 * \param c The alphabet letter that you wish to assign to.
 * \param in The `rank` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename alphabet_type>
constexpr alphabet_type & assign_rank(alphabet_type & c, underlying_rank_t<alphabet_type> const in)
    requires requires (alphabet_type c) { { c.assign_rank(uint8_t{0}) } -> alphabet_type &; }
{
    return c.assign_rank(in);
}

/*!\brief Free function wrapper that calls `.assign_rank(in)` on the first argument.
 * \tparam alphabet_type Must provide an `.assign_rank()` member function.
 * \param c An alphabet letter temporary.
 * \param in The `rank` value you wish to assign.
 * \returns The assignment result as a temporary.
 * \details
 * Use this e.g. to newly create alphabet letters from rank:
 * ~~~{.cpp}
 * auto l = assign_rank(dna5{}, 1);  // l is of type dna5 and == dna5::C
 * ~~~
 */
template <typename alphabet_type>
constexpr alphabet_type && assign_rank(alphabet_type && c, underlying_rank_t<alphabet_type> const in)
    requires requires (alphabet_type c) { { c.assign_rank(uint8_t{0}) } -> alphabet_type &; }
{
    return std::move(c.assign_rank(in));
}
//!\}
//!\}

// ============================================================================
// alphabet concept
// ============================================================================

/*!\brief The generic alphabet concept that covers most data types used in ranges.
 * \ingroup alphabet
 *
 * \details
 *
 * TODO
 */
template <typename t>
concept bool alphabet_concept = requires (t t1, t t2)
{
    // StL concepts
    requires std::is_pod_v<t> == true;
    requires std::is_swappable_v<t> == true;

    // static data members
    alphabet_size<t>::value;

    // conversion to char and rank
    { to_char(t1) } -> underlying_char_t<t>;
    { to_rank(t1) } -> underlying_rank_t<t>;
    { std::cout << t1 };

    // assignment from char and rank
    { assign_char(t1,  0) } -> t &;
    { assign_rank(t1,  0) } -> t &;
    { assign_char(t{}, 0) } -> t &&;
    { assign_rank(t{}, 0) } -> t &&;

    // required comparison operators
    { t1 == t2 } -> bool;
    { t1 != t2 } -> bool;
    { t1 <  t2 } -> bool;
    { t1 >  t2 } -> bool;
    { t1 <= t2 } -> bool;
    { t1 >= t2 } -> bool;
};

} // namespace seqan3
