// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Core alphabet concept and free function/metafunction wrappers.
 */

#pragma once

#include <iostream>

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/metafunction/function.hpp>
#include <seqan3/alphabet/adaptation/pre.hpp>
#include <seqan3/alphabet/concept_pre.hpp>
#include <seqan3/alphabet/detail/member_exposure.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{

// ------------------------------------------------------------------
// semi_alphabet_concept
// ------------------------------------------------------------------

/*!\interface seqan3::semi_alphabet_concept <>
 * \brief The basis for seqan3::alphabet_concept, but requires only rank interface (not char).
 * \ingroup alphabet
 * \extends std::Regular
 * \extends std::StrictTotallyOrdered
 *
 * This concept represents "one half" of the seqan3::alphabet_concept, it requires no
 * `char` representation and corresponding interfaces. It is mostly used internally and
 * in the composition of alphabet types (see seqan3::cartesian_composition).
 *
 * Beyond the requirements stated below, the type needs to satisfy the following standard library
 * concepts:
 *
 *   * std::Regular ("copyable and default-constructible")
 *   * std::StrictTotallyOrdered ("has all comparison operators")
 *
 * For the purpose of concept checking the types `t &` and `t &&` are also considered to satisfy
 * seqan3::semi_alphabet_concept if the type `t` satisfies it.
 *
 * It is recommended that alphabets also model seqan3::standard_layout_concept and seqan3::trivially_copyable_concept
 * and all alphabets shipped with SeqAn3 do so.
 *
 * \par Serialisation
 *
 * Types that satisfy the concept (and all refinements) can be serialised via SeqAn3
 * serialisation support.
 * The rank value is (de-)serialised, types need not provide any overloads themselves.
 *
 * \par Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT semi_alphabet_concept = std::Regular<std::remove_reference_t<t>> &&
                                std::StrictTotallyOrdered<t> &&
                                requires (t v)
{
    // static data members
    alphabet_size<t>::value;
    alphabet_size_v<t>;

    // conversion to rank
    requires           noexcept(to_rank(v));
    requires std::Same<decltype(to_rank(v)), underlying_rank_t<t>>;

    // assignment from rank
    requires           noexcept(assign_rank(v,                            0));
    requires           noexcept(assign_rank(std::remove_reference_t<t>{}, 0));
    requires std::Same<decltype(assign_rank(v,                            0)), std::remove_reference_t<t> &>;
    requires std::Same<decltype(assign_rank(std::remove_reference_t<t>{}, 0)), std::remove_reference_t<t>  >;
};
//!\endcond

// ------------------------------------------------------------------
// alphabet_concept
// ------------------------------------------------------------------

/*!\interface seqan3::alphabet_concept <>
 * \extends seqan3::semi_alphabet_concept
 * \brief The generic alphabet concept that covers most data types used in ranges.
 * \ingroup alphabet
 *
 * This is the core alphabet concept that many other alphabet concepts refine.
 *
 * It defines the requirements for the rank interface and the character interface,
 * as well as the requirement for size and comparability. For more details, see
 * the \ref alphabet module.
 *
 * For the purpose of concept checking the types `t &` and `t &&` are also considered to satisfy
 * seqan3::alphabet_concept if the type `t` satisfies it.
 *
 * \par Serialisation
 *
 * Types that satisfy the concept (and all refinements) can be serialised via SeqAn3
 * serialisation support.
 * The rank value is (de-)serialised, types need not provide any overloads themselves.
 *
 * \par Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT alphabet_concept = semi_alphabet_concept<t> && requires (t v)
{
    // conversion to char
    requires           noexcept(to_char(v));
    requires std::Same<decltype(to_char(v)), underlying_char_t<t>>;

    // assignment from char
    requires           noexcept(assign_char(v,                            0));
    requires           noexcept(assign_char(std::remove_reference_t<t>{}, 0));
    requires std::Same<decltype(assign_char(v,                            0)), std::remove_reference_t<t> &>;
    requires std::Same<decltype(assign_char(std::remove_reference_t<t>{}, 0)), std::remove_reference_t<t>  >;

    // chars can be checked for validity
    requires std::Same<decltype(char_is_valid_for<t>(char{0})), bool>;

    // strict assignment from char
    requires std::Same<decltype(assign_char_strict(v,                            0)), std::remove_reference_t<t> &>;
    requires std::Same<decltype(assign_char_strict(std::remove_reference_t<t>{}, 0)), std::remove_reference_t<t>  >;
};
//!\endcond

// ------------------------------------------------------------------
//  serialisation
// ------------------------------------------------------------------

/*!\cond DEV
 * \name Generic serialisation functions for all seqan3::semi_alphabet_concept
 * \brief All types that satisfy seqan3::semi_alphabet_concept can be serialised via Cereal.
 *
 * \{
 */
/*!
 * \brief Save an alphabet letter to stream.
 * \tparam archive_t Must satisfy seqan3::cereal_output_archive_concept.
 * \tparam alphabet_t Type of l; must satisfy seqan3::semi_alphabet_concept.
 * \param l The alphabet letter.
 * \relates seqan3::semi_alphabet_concept
 *
 * \details
 *
 * Delegates to seqan3::semi_alphabet_concept::to_rank.
 *
 * \attention These functions are never called directly, see the \ref alphabet module on how to use serialisation.
 */
template <cereal_output_archive_concept archive_t, semi_alphabet_concept alphabet_t>
underlying_rank_t<alphabet_t> CEREAL_SAVE_MINIMAL_FUNCTION_NAME(archive_t const &, alphabet_t const & l)
{
    return to_rank(l);
}

/*!\brief Restore an alphabet letter from a saved rank.
 * \tparam archive_t Must satisfy seqan3::cereal_input_archive_concept.
 * \tparam wrapped_alphabet_t A seqan3::semi_alphabet_concept after Cereal mangles it up.
 * \param l The alphabet letter (cereal wrapped).
 * \param r The assigned value.
 * \relates seqan3::semi_alphabet_concept
 *
 * \details
 *
 * Delegates to seqan3::semi_alphabet_concept::assign_rank.
 *
 * \attention These functions are never called directly, see the \ref alphabet module on how to use serialisation.
 */
template <cereal_input_archive_concept archive_t, typename wrapped_alphabet_t>
void CEREAL_LOAD_MINIMAL_FUNCTION_NAME(archive_t const &,
                                       wrapped_alphabet_t && l,
                                       underlying_rank_t<detail::strip_cereal_wrapper_t<wrapped_alphabet_t>> const & r)
    requires semi_alphabet_concept<detail::strip_cereal_wrapper_t<wrapped_alphabet_t>>
{
    assign_rank(static_cast<detail::strip_cereal_wrapper_t<wrapped_alphabet_t> &>(l), r);
}
/*!\}
 * \endcond
 */

} // namespace seqan3

namespace seqan3::detail
{
// ------------------------------------------------------------------
// constexpr_semi_alphabet_concept
// ------------------------------------------------------------------

/*!\interface seqan3::detail::constexpr_semi_alphabet_concept <>
 * \brief A seqan3::semi_alphabet_concept that has a constexpr default constructor and constexpr accessors.
 * \ingroup alphabet
 * \extends seqan3::semi_alphabet_concept
 *
 * The same as seqan3::semi_alphabet_concept, except that all required functions are also required to be callable
 * in a `constexpr`-context.
 *
 * \par Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT constexpr_semi_alphabet_concept = semi_alphabet_concept<t> && requires
{
    // currently only tests rvalue interfaces, because we have no constexpr values in this scope to get references to
    requires SEQAN3_IS_CONSTEXPR(to_rank(std::remove_reference_t<t>{}));
    requires SEQAN3_IS_CONSTEXPR(assign_rank(std::remove_reference_t<t>{}, 0));
};
//!\endcond

// ------------------------------------------------------------------
// constexpr_alphabet_concept
// ------------------------------------------------------------------

/*!\interface seqan3::detail::constexpr_alphabet_concept <>
 * \brief A seqan3::alphabet_concept that has constexpr accessors.
 * \ingroup alphabet
 * \extends seqan3::detail::constexpr_semi_alphabet_concept
 * \extends seqan3::alphabet_concept
 *
 * The same as seqan3::alphabet_concept, except that the following interface requirements are also required to be
 * callable in a `constexpr`-context:
 *
 *   * seqan3::alphabet_concept::to_char
 *   * seqan3::alphabet_concept::assign_char
 *   * seqan3::alphabet_concept::char_is_valid_for
 *
 * The only exception is:
 *
 *   * seqan3::alphabet_concept::assign_char_strict
 *
 * \par Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT constexpr_alphabet_concept = constexpr_semi_alphabet_concept<t> &&
                                     alphabet_concept<t> &&
                                     requires
{
    // currently only tests rvalue interfaces, because we have no constexpr values in this scope to get references to
    requires SEQAN3_IS_CONSTEXPR(to_char(std::remove_reference_t<t>{}));
    requires SEQAN3_IS_CONSTEXPR(assign_char(std::remove_reference_t<t>{}, underlying_char_t<t>{}));
    requires SEQAN3_IS_CONSTEXPR(char_is_valid_for<t>(char{0}));
};
//!\endcond

} // namespace seqan3::detail

#include <seqan3/alphabet/detail/hash.hpp>
