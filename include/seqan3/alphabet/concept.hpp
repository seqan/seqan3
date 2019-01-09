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
//TODO(rrahn): Change to template <typename t2, typename t = std::remove_reference_t<t2>>
// in order to get rid of the remove_reference_t within the concept, after the ICE
// get's fixed. See issue #228
template <typename t>
concept semi_alphabet_concept = std::Regular<std::remove_reference_t<t>> &&
                                std::StrictTotallyOrdered<t> &&
                                requires (t t1, t t2)
{
    // static data members
    alphabet_size<std::remove_reference_t<t>>::value;
    alphabet_size_v<std::remove_reference_t<t>>;

    // conversion to rank
    requires noexcept(to_rank(t1));
    requires std::Same<decltype(to_rank(t1)), underlying_rank_t<std::remove_reference_t<t>>>;

    // assignment from rank
    requires noexcept(assign_rank(t1,  0));
    requires std::Same<decltype(assign_rank(t1,  0)), std::remove_reference_t<t> &>;
    requires noexcept(assign_rank(std::remove_reference_t<t>{}, 0));
    requires std::Same<decltype(assign_rank(std::remove_reference_t<t>{}, 0)), std::remove_reference_t<t> &&>;
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
//TODO(rrahn): Change to template <typename t2, typename t = std::remove_reference_t<t2>>
// in order to get rid of the remove_reference_t within the concept, after the ICE
// get's fixed. See issue #228
template <typename t>
concept alphabet_concept = semi_alphabet_concept<t> && requires (t t1, t t2)
{
    // conversion to char
    requires noexcept(to_char(t1));
    requires std::Same<decltype(to_char(t1)), underlying_char_t<std::remove_reference_t<t>>>;

    // assignment from char
    requires noexcept(assign_char(t1,  0));
    requires std::Same<decltype(assign_char(t1,  0)), std::remove_reference_t<t> &>;
    requires noexcept(assign_char(std::remove_reference_t<t>{},  0));
    requires std::Same<decltype(assign_char(std::remove_reference_t<t>{}, 0)), std::remove_reference_t<t> &&>;
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
    assign_rank(static_cast<detail::strip_cereal_wrapper_t<wrapped_alphabet_t>&&>(l), r);
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
 * The same as seqan3::semi_alphabet_concept, except that all required functions are also required to be
 * `constexpr`-qualified.
 *
 * \par Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
concept constexpr_semi_alphabet_concept = semi_alphabet_concept<t> && requires
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
 * The same as seqan3::alphabet_concept, except that all required functions are also required to be
 * `constexpr`-qualified.
 *
 * \par Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
concept constexpr_alphabet_concept = constexpr_semi_alphabet_concept<t> &&
                                     alphabet_concept<t> &&
                                     requires
{
    // currently only tests rvalue interfaces, because we have no constexpr values in this scope to get references to
    requires SEQAN3_IS_CONSTEXPR(to_char(std::remove_reference_t<t>{}));
    requires SEQAN3_IS_CONSTEXPR(assign_char(std::remove_reference_t<t>{},
                                             underlying_char_t<std::remove_reference_t<t>>{}));
};
//!\endcond

} // namespace seqan3::detail

#include <seqan3/alphabet/detail/hash.hpp>
