// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Provides seqan3::AminoacidAlphabet.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>

// ============================================================================
// concept
// ============================================================================

namespace seqan3
{

/*!\brief Identifies amino acid alphabets.
 * \implements seqan3::UnaryTypeTrait
 * \ingroup aminoacid
 *
 * \details
 *
 * Since an amino acid alphabet has no specific characteristics (like the complement
 * function for nucleotide alphabets), we distinguish an amino acid alphabet by
 * the seqan3::is_aminoacid unary type trait.
 *
 * ### Customisation point
 *
 * If you define your own alphabet and want it to be recognised as an amino acid
 * alphabet by SeqAn, you need to specialise this type trait for your type and
 * have it inherit std::true_type.
 *
 * \include test/snippet/alphabet/aminoacid/is_aminoacid.cpp
 */
template <typename type>
struct is_aminoacid : std::false_type {};

//!\brief Helper variable that delegates to seqan3::is_aminoacid<type>::value (UnaryTypeTrait shortcut).
//!\relates seqan3::is_aminoacid
//!\ingroup aminoacid
template <typename type>
constexpr bool is_aminoacid_v = is_aminoacid<type>::value;

/*!\interface seqan3::AminoacidAlphabet <>
 * \extends seqan3::Alphabet
 * \brief A concept that indicates whether an alphabet represents amino acids.
 * \ingroup aminoacid
 *
 * Since an amino acid alphabet has no specific characteristics (like the complement
 * function for nucleotide alphabets), we distinguish an amino acid alphabet by
 * the seqan3::is_aminoacid type trait.
 *
 * ###Concepts and doxygen
 * The requirements for this concept are given as related functions and type traits.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename type>
SEQAN3_CONCEPT AminoacidAlphabet = Alphabet<type> && is_aminoacid_v<remove_cvref_t<type>>;
//!\endcond

} // namespace seqan3
