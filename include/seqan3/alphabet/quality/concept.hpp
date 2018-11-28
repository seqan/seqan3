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
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \brief Quality alphabet concept.
 */

#pragma once

#include <iostream>
#include <string>

#include <seqan3/alphabet/concept.hpp>

namespace seqan3
{

// ------------------------------------------------------------------
// Member exposure for the seqan3::quality_concept
// ------------------------------------------------------------------

/*!\name Helpers for seqan3::quality_concept
 * \brief These functions and metafunctions expose member variables and types so
 * that the type can model the seqan3::quality_concept.
 * \ingroup quality
 * \{
 */

template <typename alphabet_type>
struct underlying_phred
{};

/*!\brief The internal phred type.
 * \ingroup quality
 * \tparam alphabet_type The type of alphabet. Must model the seqan3::quality_concept.
 *
 * The underlying_phred type requires the quality_concept.
 */
template <typename alphabet_with_member_type>
//!\cond
    requires requires () { typename alphabet_with_member_type::phred_type; }
//!\endcond
struct underlying_phred<alphabet_with_member_type>
{
    //!\brief The underlying phred data type.
    using type = typename alphabet_with_member_type::phred_type;
};

/*!\brief The internal phred type.
 * \ingroup quality
 * \tparam alphabet_type The type of alphabet. Must model the seqan3::quality_concept.
 *
 * The underlying_phred type requires the quality_concept.
 */
template <typename alphabet_type>
using underlying_phred_t = typename underlying_phred<alphabet_type>::type;

/*!\brief The public setter function of a phred score.
 * \ingroup quality
 * \tparam    alphabet_type The type of alphabet. Must model the seqan3::quality_concept.
 * \param[in] chr           The quality value to assign a score.
 * \param[in] in            The character to representing the phred score.
 *
 * The underlying_phred type requires the quality_concept.
 */
template <typename alphabet_type>
//!\cond
    requires requires (alphabet_type v) { { v.assign_phred('c') }; }
//!\endcond
constexpr alphabet_type & assign_phred(alphabet_type & chr, char const in)
{
    return chr.assign_phred(in);
}

template <typename alphabet_type>
//!\cond
    requires requires (alphabet_type v) { { v.assign_phred('c') }; }
//!\endcond
constexpr alphabet_type assign_phred(alphabet_type && chr, char const in)
{
    return chr.assign_phred(in);
}

/*!\brief The public getter function for the phred representation of a score.
 * \ingroup quality
 * \tparam    alphabet_type The type of alphabet. Must model the seqan3::quality_concept.
 * \param[in] chr           The quality value to convert into the phred score.
 *
 * The underlying_phred type requires the quality_concept.
 */
template <typename alphabet_type>
//!\cond
    requires requires (alphabet_type v) { { v.to_phred() }; }
//!\endcond
constexpr underlying_phred_t<alphabet_type> to_phred(alphabet_type const & chr)
{
    return chr.to_phred();
}
//\}

// ------------------------------------------------------------------
// seqan3::quality_concept
// ------------------------------------------------------------------

/*!\interface seqan3::quality_concept <>
 * \extends seqan3::alphabet_concept
 * \brief A concept that indicates whether an alphabet represents quality scores.
 * \ingroup quality
 *
 * In addition to the requirements for seqan3::alphabet_concept, the
 * quality_concept introduces a requirement for conversion functions from and to
 * a Phred score.
 *
 * \par Concepts and doxygen
 * The requirements for this concept are given as related functions and
 * metafunctions. Types that satisfy this concept are shown as "implementing
 * this interface".
 */
//!\cond
template<typename q>
concept quality_concept = requires(q quality)
{
    requires alphabet_concept<q>;

    { assign_phred(quality, typename q::rank_type{}) } -> q;
    { to_phred(quality) } -> const typename q::phred_type;
    typename underlying_phred<q>::type;
};
//!\endcond

} // namespace seqan
