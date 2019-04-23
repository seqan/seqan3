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
// Semialphabet
// ------------------------------------------------------------------

/*!\interface seqan3::Semialphabet <>
 * \brief The basis for seqan3::Alphabet, but requires only rank interface (not char).
 * \ingroup alphabet
 * \extends std::Regular
 * \extends std::StrictTotallyOrdered
 *
 * This concept represents "one half" of the seqan3::Alphabet, it requires no
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
 * seqan3::Semialphabet if the type `t` satisfies it.
 *
 * It is recommended that alphabets also model seqan3::StandardLayout and seqan3::TriviallyCopyable
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
SEQAN3_CONCEPT Semialphabet = std::Regular<std::remove_reference_t<t>> &&
                                std::StrictTotallyOrdered<t> &&
                                requires (t v)
{
    // static data members
    alphabet_size<t>::value;
    alphabet_size_v<t>;

    // conversion to rank
    requires           noexcept(to_rank(v));
    requires std::Same<decltype(to_rank(v)), alphabet_rank_t<t>>;

    // assignment from rank
    requires           noexcept(assign_rank_to(0, v                           ));
    requires           noexcept(assign_rank_to(0, std::remove_reference_t<t>{}));
    requires std::Same<decltype(assign_rank_to(0, v                           )), std::remove_reference_t<t> &>;
    requires std::Same<decltype(assign_rank_to(0, std::remove_reference_t<t>{})), std::remove_reference_t<t>  >;
};
//!\endcond

// ------------------------------------------------------------------
// Alphabet
// ------------------------------------------------------------------

/*!\interface seqan3::Alphabet <>
 * \extends seqan3::Semialphabet
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
 * seqan3::Alphabet if the type `t` satisfies it.
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
SEQAN3_CONCEPT Alphabet = Semialphabet<t> && requires (t v)
{
    // conversion to char
    requires           noexcept(to_char(v));
    requires std::Same<decltype(to_char(v)), alphabet_char_t<t>>;

    // assignment from char
    requires           noexcept(assign_char_to(0, v                           ));
    requires           noexcept(assign_char_to(0, std::remove_reference_t<t>{}));
    requires std::Same<decltype(assign_char_to(0, v                           )), std::remove_reference_t<t> &>;
    requires std::Same<decltype(assign_char_to(0, std::remove_reference_t<t>{})), std::remove_reference_t<t>  >;

    // chars can be checked for validity
    requires std::Same<decltype(char_is_valid_for<t>(char{0})), bool>;

    // strict assignment from char
    requires std::Same<decltype(assign_char_strictly_to(0, v                           )), std::remove_reference_t<t> &>;
    requires std::Same<decltype(assign_char_strictly_to(0, std::remove_reference_t<t>{})), std::remove_reference_t<t>  >;
};
//!\endcond

// ------------------------------------------------------------------
//  serialisation
// ------------------------------------------------------------------

/*!\cond DEV
 * \name Generic serialisation functions for all seqan3::Semialphabet
 * \brief All types that satisfy seqan3::Semialphabet can be serialised via Cereal.
 *
 * \{
 */
/*!
 * \brief Save an alphabet letter to stream.
 * \tparam archive_t Must satisfy seqan3::CerealOutputArchive.
 * \tparam alphabet_t Type of l; must satisfy seqan3::Semialphabet.
 * \param l The alphabet letter.
 * \relates seqan3::Semialphabet
 *
 * \details
 *
 * Delegates to seqan3::Semialphabet::to_rank.
 *
 * \attention These functions are never called directly, see the \ref alphabet module on how to use serialisation.
 */
template <CerealOutputArchive archive_t, Semialphabet alphabet_t>
alphabet_rank_t<alphabet_t> CEREAL_SAVE_MINIMAL_FUNCTION_NAME(archive_t const &, alphabet_t const & l)
{
    return to_rank(l);
}

/*!\brief Restore an alphabet letter from a saved rank.
 * \tparam archive_t Must satisfy seqan3::CerealInputArchive.
 * \tparam wrapped_alphabet_t A seqan3::Semialphabet after Cereal mangles it up.
 * \param l The alphabet letter (cereal wrapped).
 * \param r The assigned value.
 * \relates seqan3::Semialphabet
 *
 * \details
 *
 * Delegates to seqan3::Semialphabet::assign_rank.
 *
 * \attention These functions are never called directly, see the \ref alphabet module on how to use serialisation.
 */
template <CerealInputArchive archive_t, typename wrapped_alphabet_t>
void CEREAL_LOAD_MINIMAL_FUNCTION_NAME(archive_t const &,
                                       wrapped_alphabet_t && l,
                                       alphabet_rank_t<detail::strip_cereal_wrapper_t<wrapped_alphabet_t>> const & r)
    requires Semialphabet<detail::strip_cereal_wrapper_t<wrapped_alphabet_t>>
{
    assign_rank_to(r, static_cast<detail::strip_cereal_wrapper_t<wrapped_alphabet_t> &>(l));
}
/*!\}
 * \endcond
 */

} // namespace seqan3

namespace seqan3::detail
{
// ------------------------------------------------------------------
// ConstexprSemialphabet
// ------------------------------------------------------------------

/*!\interface seqan3::detail::ConstexprSemialphabet <>
 * \brief A seqan3::Semialphabet that has a constexpr default constructor and constexpr accessors.
 * \ingroup alphabet
 * \extends seqan3::Semialphabet
 *
 * The same as seqan3::Semialphabet, except that all required functions are also required to be callable
 * in a `constexpr`-context.
 *
 * \par Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT ConstexprSemialphabet = Semialphabet<t> && requires
{
    // currently only tests rvalue interfaces, because we have no constexpr values in this scope to get references to
    requires SEQAN3_IS_CONSTEXPR(to_rank(std::remove_reference_t<t>{}));
    requires SEQAN3_IS_CONSTEXPR(assign_rank_to(0, std::remove_reference_t<t>{}));
};
//!\endcond

// ------------------------------------------------------------------
// ConstexprAlphabet
// ------------------------------------------------------------------

/*!\interface seqan3::detail::ConstexprAlphabet <>
 * \brief A seqan3::Alphabet that has constexpr accessors.
 * \ingroup alphabet
 * \extends seqan3::detail::ConstexprSemialphabet
 * \extends seqan3::Alphabet
 *
 * The same as seqan3::Alphabet, except that the following interface requirements are also required to be
 * callable in a `constexpr`-context:
 *
 *   * seqan3::Alphabet::to_char
 *   * seqan3::Alphabet::assign_char
 *   * seqan3::Alphabet::char_is_valid_for
 *
 * The only exception is:
 *
 *   * seqan3::Alphabet::assign_char_strict
 *
 * \par Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT ConstexprAlphabet = ConstexprSemialphabet<t> &&
                                     Alphabet<t> &&
                                     requires
{
    // currently only tests rvalue interfaces, because we have no constexpr values in this scope to get references to
    requires SEQAN3_IS_CONSTEXPR(to_char(std::remove_reference_t<t>{}));
    requires SEQAN3_IS_CONSTEXPR(assign_char_to(alphabet_char_t<t>{}, std::remove_reference_t<t>{}));
    requires SEQAN3_IS_CONSTEXPR(char_is_valid_for<t>(char{0}));
};
//!\endcond

} // namespace seqan3::detail

#include <seqan3/alphabet/detail/hash.hpp>
