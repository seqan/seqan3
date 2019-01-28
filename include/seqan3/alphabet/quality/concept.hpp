// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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
SEQAN3_CONCEPT quality_concept = requires(q quality)
{
    requires alphabet_concept<q>;

    { assign_phred(quality, typename q::rank_type{}) } -> q;
    { to_phred(quality) } -> const typename q::phred_type;
    typename underlying_phred<q>::type;
};
//!\endcond

} // namespace seqan
