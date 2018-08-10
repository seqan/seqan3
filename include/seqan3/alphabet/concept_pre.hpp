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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief seqan3::alphabet_concept metafunction base classes.
 *
 * Include this file, if you implement an alphabet type with free/global function
 * and metafunction interfaces.
 *
 * \attention
 *
 * Note that you need to strictly follow this include order:
 * ```cpp
 * #include <alphabet/concept_pre.hpp>
 *
 * // your custom alphabet
 *
 * #include <alphabet/concept.hpp>
 * ```
 *
 * If you include `concept.hpp` before your definitions, than your type will
 * not be resolved as satisfying seqan3::alphabet_concept.
 *
 * This is not true for custom alphabets implementing the interfaces as
 * member functions/variables/types.
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{

// ------------------------------------------------------------------
// seqan3::semi_alphabet_concept
// ------------------------------------------------------------------

/*!\name Requirements for seqan3::semi_alphabet_concept
 * \brief You can expect these functions on all types that implement seqan3::semi_alphabet_concept.
 * \relates seqan3::semi_alphabet_concept
 * \{
 */
/*!\brief The `rank_type` of the semi_alphabet. [type metafunction base template]
 * \tparam semi_alphabet_type Must satisfy seqan3::semi_alphabet_concept.
 * \ingroup alphabet
 *
 * \par Helper template alias
 *   seqan3::underlying_rank_t as a shorthand for `typename seqan3::underlying_rank<semi_alphabet_type>::%type`
 *
 * \attention This is the base template, it needs to be specialised.
 */
template <typename semi_alphabet_type>
struct underlying_rank{};

/*!\brief The `rank_type` of the semi_alphabet. [type metafunction shortcut]
 * \ingroup alphabet
 *
 * \attention Do not specialise this shortcut, instead specialise seqan3::underlying_rank.
 */
template <typename semi_alphabet_type>
using underlying_rank_t = typename underlying_rank<semi_alphabet_type>::type;

/*!\brief The size of the alphabet. [value metafunction base template]
 * \tparam alphabet_type Must satisfy seqan3::semi_alphabet_concept.
 * \ingroup alphabet
 *
 * This is the expression to retrieve the value:
 * ```cpp
 * auto i = seqan3::alphabet_size<alphabet_type>::value;
 * // or
 * auto i = seqan3::alphabet_size_v<alphabet_type>;
 * ```
 * The type of the variable is seqan3::underlying_rank_t<alphabet_type>.
 *
 * \par Helper variable template
 *   seqan3::alphabet_size_v as a shorthand for `seqan3::alphabet_size<alphabet_type>::%value`
 *
 * \attention This is the base template, it needs to be specialised.
 */
template <typename alphabet_type>
struct alphabet_size{};

/*!\brief The size of the alphabet. [value metafunction shortcut]
 * \tparam alphabet_type Must satisfy seqan3::semi_alphabet_concept.
 * \ingroup alphabet
 *
 * \attention Do not specialise this shortcut, instead specialise seqan3::alphabet_size.
 */
template <typename alphabet_type>
//!\cond
    requires requires (alphabet_type alph) { alphabet_size<alphabet_type>::value; }
//!\endcond
constexpr auto alphabet_size_v = alphabet_size<alphabet_type>::value;

/*!\fn rank_type seqan3::to_rank(semi_alphabet_concept const alph)
 * \brief Returns the alphabet letter's value in rank representation.
 * \ingroup alphabet
 * \param alph The alphabet letter that you wish to convert to rank.
 * \returns The letter's value in the alphabet's rank type (seqan3::underlying_rank).
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
// just implement the interface

/*!\fn semi_alphabet_concept && seqan3::assign_rank(semi_alphabet_concept && alph, rank_type const rank)
 * \brief Returns the alphabet letter's value in rank representation.
 * \ingroup alphabet
 * \param alph The alphabet letter that you wish to assign to.
 * \param rank The rank you wish to assign.
 * \returns A reference to `alph` or a temporary if `alph` was a temporary.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
// just implement the interface
//!\}

// ------------------------------------------------------------------
// seqan3::alphabet_concept
// ------------------------------------------------------------------

/*!\name Requirements for seqan3::alphabet_concept
 * \brief You can expect these functions on all types that implement seqan3::alphabet_concept.
 * \relates seqan3::alphabet_concept
 * \{
 */

/*!\brief The `char_type` of the alphabet. [type metafunction base template]
 * \tparam alphabet_type Must satisfy seqan3::alphabet_concept.
 * \ingroup alphabet
 *
 * \par Helper template alias
 *   seqan3::underlying_char_t as a shorthand for `typename seqan3::underlying_char<alphabet_type>::%type`
 *
 * \attention This is the base template, it needs to be specialised.
 */
template <typename alphabet_type>
struct underlying_char{};

/*!\brief The `char_type` of the alphabet. [type metafunction shortcut]
 * \ingroup alphabet
 *
 * \attention Do not specialise this shortcut, instead specialise seqan3::underlying_char.
 */
template <typename alphabet_type>
using underlying_char_t = typename underlying_char<alphabet_type>::type;

/*!\fn char_type seqan3::to_char(alphabet_concept const alph)
 * \brief Returns the alphabet letter's value in character representation.
 * \ingroup alphabet
 * \param alph The alphabet letter that you wish to convert to char.
 * \returns The letter's value in the alphabet's char type (seqan3::underlying_char).
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
// just implement the interface

/*!\fn alphabet_concept && seqan3::assign_char(alphabet_concept && alph, char_type const chr)
 * \brief Returns the alphabet letter's value in character representation.
 * \ingroup alphabet
 * \param alph The alphabet letter that you wish to assign to.
 * \param chr The `char` you wish to assign.
 * \returns A reference to `alph` or a temporary if `alph` was a temporary.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
// just implement the interface

/*!\fn std::ostream & operator<<(std::ostream & os, alphabet_type const alph)
 * \brief Ostream operator for the alphabet.
 * \ingroup alphabet
 * \param os The output stream you are printing to.
 * \param alph The alphabet letter that you wish to convert to char.
 * \returns A reference to the output stream.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
// just implement the interface

//!\}

// ------------------------------------------------------------------
// seqan3::rna_structure_concept
// ------------------------------------------------------------------

/*!\name Requirements for seqan3::rna_structure_concept
 * \brief You can expect these functions on all types that implement seqan3::rna_structure_concept.
 * \relates seqan3::rna_structure_concept
 * \{
 */
/*!\brief Metafunction that indicates to what extent an alphabet can handle pseudoknots.
 * [value metafunction base template]
 * \tparam alphabet_type The alphabet type whose pseudoknot ability is queried.
 * \ingroup alphabet
 * \details The value is the maximum allowed depth of pseudoknots.
 * A value of 1 denotes no pseudoknots `((....))`,
 * while higher values denote the maximum allowed complexity of
 * crossing interactions, e.g. depth 2 `(({....))}` or depth 3 `({[....)}]`.
 *
 * This is the expression to retrieve the value:
 * ```cpp
 * uint8_t pk_support = seqan3::pseudoknot_support<alphabet_type>::value;
 * // or
 * uint8_t pk_support = seqan3::pseudoknot_support_v<alphabet_type>;
 * ```
 *
 * \par Helper variable template
 *   seqan3::pseudoknot_support_v as a shorthand for `seqan3::pseudoknot_support<alphabet_type>::%value`
 *
 * \attention This is the base template, it needs to be specialised.
 */
template<typename alphabet_type>
struct pseudoknot_support{};

/*!\brief The pseudoknot ability of the alphabet. [value metafunction shortcut]
 * \ingroup alphabet
 *
 * \attention Do not specialise this shortcut, instead specialise seqan3::pseudoknot_support.
 */
template<typename alphabet_type>
constexpr uint8_t pseudoknot_support_v = pseudoknot_support<alphabet_type>::value;
//!\}

} // namespace seqan3

// ------------------------------------------------------------------
// seqan3::char_adaption_concept
// ------------------------------------------------------------------

namespace seqan3::detail
{

/*!\brief Metafunction that indicates whether a type is a char alphabet adaptation.
 * \ingroup alphabet
 *
 * \par Helper variable template
 *   seqan3::is_char_adaptation_v as a shorthand for `seqan3::is_char_adaptation<type>::%value`
 */
template <typename type>
struct is_char_adaptation :
    std::false_type
{};

//!\brief Shortcut for seqan3::detail::is_char_adaptation.
template <typename type>
constexpr bool is_char_adaptation_v = is_char_adaptation<type>::value;

// ------------------------------------------------------------------
// seqan3::uint_adaption_concept
// ------------------------------------------------------------------

//!\brief Metafunction that indicates whether a type is a uint alphabet adaptation.
template <typename type>
struct is_uint_adaptation :
    std::false_type
{};

//!\brief Shortcut for seqan3::detail::is_uint_adaptation.
template <typename type>
constexpr bool is_uint_adaptation_v = is_uint_adaptation<type>::value;

} // namespace seqan3::detail
