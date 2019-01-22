// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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
 * \snippet test/snippet/alphabet/concept_pre.cpp include order
 *
 * If you include `concept.hpp` before your definitions, than your type will
 * not be resolved as satisfying seqan3::alphabet_concept.
 *
 * This is not `true` for custom alphabets implementing the interfaces as
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
struct underlying_rank
{};

template <typename semi_alphabet_type>
struct underlying_rank<semi_alphabet_type &> : underlying_rank<semi_alphabet_type>
{};

template <typename semi_alphabet_type>
struct underlying_rank<semi_alphabet_type &&> : underlying_rank<semi_alphabet_type>
{};

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
 * \snippet test/snippet/alphabet/concept_pre.cpp value retrieval
 * The type of the variable is seqan3::underlying_rank_t<alphabet_type>.
 *
 * \par Helper variable template
 *   seqan3::alphabet_size_v as a shorthand for `seqan3::alphabet_size<alphabet_type>::%value`
 *
 * \attention This is the base template, it needs to be specialised.
 */
template <typename alphabet_type>
struct alphabet_size
{};

template <typename alphabet_type>
struct alphabet_size<alphabet_type &> : alphabet_size<alphabet_type>
{};

template <typename alphabet_type>
struct alphabet_size<alphabet_type &&> : alphabet_size<alphabet_type>
{};

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

/*!\fn rank_type seqan3::to_rank(semi_alphabet_type const alph) noexcept
 * \brief Returns the alphabet letter's value in rank representation.
 * \ingroup alphabet
 * \param alph The alphabet letter that you wish to convert to rank.
 * \returns The letter's value in the alphabet's rank type (seqan3::underlying_rank).
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 *
 * This function is required to not throw and be marked as `noexcept`.
 */
// just implement the interface

/*!\fn semi_alphabet_type && seqan3::assign_rank(semi_alphabet_type && alph, rank_type const rank) noexcept
 * \brief Returns the alphabet letter's value in rank representation.
 * \ingroup alphabet
 * \param alph The alphabet letter that you wish to assign to.
 * \param rank The rank you wish to assign.
 * \returns A reference to `alph` or a temporary if `alph` was a temporary.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 *
 *   * This function is required to not throw and be marked as `noexcept`.
 *   * Assigning rank values that are larger than or equal to the alphabet's size may be undefined behaviour.
 *   * Implementations that wish to avoid undefined behaviour cannot throw exceptions, but may use contracts or
 *     assertions, or simply convert all invalid rank values to valid ones.
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
struct underlying_char
{};

template <typename alphabet_type>
struct underlying_char<alphabet_type &> : underlying_char<alphabet_type>
{};

template <typename alphabet_type>
struct underlying_char<alphabet_type &&> : underlying_char<alphabet_type>
{};

/*!\brief The `char_type` of the alphabet. [type metafunction shortcut]
 * \ingroup alphabet
 *
 * \attention Do not specialise this shortcut, instead specialise seqan3::underlying_char.
 */
template <typename alphabet_type>
using underlying_char_t = typename underlying_char<alphabet_type>::type;

/*!\fn char_type seqan3::to_char(alphabet_type const alph) noexcept
 * \brief Returns the alphabet letter's value in character representation.
 * \ingroup alphabet
 * \param alph The alphabet letter that you wish to convert to char.
 * \returns The letter's value in the alphabet's char type (seqan3::underlying_char).
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 *
 * This function is required to not throw and be marked as `noexcept`.
 */
// just implement the interface

/*!\fn alphabet_type && seqan3::assign_char(alphabet_type && alph, char_type const chr) noexcept
 * \brief Assigns a character value to an alphabet object.
 * \ingroup alphabet
 * \param alph The alphabet letter that you wish to assign to.
 * \param chr The `char` you wish to assign.
 * \returns A reference to `alph` or a temporary if `alph` was a temporary.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 *
 *   * This function is required to not throw and be marked as `noexcept`.
 *   * This function is required to accept all values of the character type (in addition to not throwing, it may
 *     not assert or prevent invalid character values via contracts).
 *   * All character values must leave the alphabet object in a valid and defined state; implementations typically
 *     convert invalid characters to valid characters before assignment.
 */
// just implement the interface

/*!\fn alphabet_type && seqan3::assign_char_strict(alphabet_type && alph, char_type const chr)
 * \brief Assigns a character value to an alphabet object.
 * \ingroup alphabet
 * \param alph The alphabet letter that you wish to assign to.
 * \param chr The `char` you wish to assign.
 * \returns A reference to `alph` or a temporary if `alph` was a temporary.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 *
 * Typical implementations of this function throw seqan3::invalid_char_assignment if seqan3::char_is_valid_for returns
 * `false` and call seqan3::assign_char otherwise, but as it is not required that seqan3::char_is_valid_for ever return
 * `false` for a given alphabet, this function may also be implemented as simply forwarding to seqan3::assign_char.
 */
// just implement the interface

/*!\fn bool seqan3::char_is_valid_for<alphabet_type>(char_type const chr) noexcept
 * \brief Whether the given character value represents a value of the given alphabet.
 * \ingroup alphabet
 * \tparam alphabet_type The alphabet type to query.
 * \param chr The character value to query.
 * \returns `true` if the given character has a one-to-one mapping to an alphabet value, `false` otherwise.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 *
 * This function is required to not throw and be marked as `noexcept`.
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
 * \ingroup structure
 * \details The value is the maximum allowed depth of pseudoknots.
 * A value of 1 denotes no pseudoknots `((....))`,
 * while higher values denote the maximum allowed complexity of
 * crossing interactions, e.g. depth 2 `(({....))}` or depth 3 `({[....)}]`.
 *
 * This is the expression to retrieve the value:
 * \snippet test/snippet/alphabet/concept_pre.cpp pseudoknot value retrieval
 *
 * \par Helper variable template
 *   seqan3::max_pseudoknot_depth_v as a shorthand for `seqan3::max_pseudoknot_depth<alphabet_type>::%value`
 *
 * \attention This is the base template, it needs to be specialised.
 */
template<typename alphabet_type>
struct max_pseudoknot_depth{};

/*!\brief The pseudoknot ability of the alphabet. [value metafunction shortcut]
 * \ingroup structure
 *
 * \attention Do not specialise this shortcut, instead specialise seqan3::max_pseudoknot_depth.
 */
template<typename alphabet_type>
//!\cond
    requires requires (alphabet_type alph) { max_pseudoknot_depth<alphabet_type>::value; }
//!\endcond
constexpr uint8_t max_pseudoknot_depth_v = max_pseudoknot_depth<alphabet_type>::value;
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
