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

/*!\file
 * \ingroup alphabet
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Core alphabet concept and free function/metafunction wrappers.
 */

#pragma once

#define SEQAN3_ALPHABET_CONCEPT_INCLUDED

#include <iostream>

#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/alphabet/concept_pre.hpp>
#include <seqan3/alphabet/detail/member_exposure.hpp>

namespace seqan3
{

/*!\interface seqan3::semi_alphabet_concept <>
 * \brief The basis for seqan3::alphabet_concept, but requires only rank interface (not char).
 * \ingroup alphabet
 *
 * This concept represents "one half" of the seqan3::alphabet_concept, it requires no
 * `char` representation and corresponding interfaces. It is mostly used internally and
 * in the composition of alphabet types (see seqan3::cartesian_composition).
 *
 * Beyond the requirements stated below, the type needs to be a plain old datatype (`std::is_pod_v`)
 * and be swappable (`std::is_swappable_v`).
 *
 * \todo add comparison operators
 *
 * \par Concepts and doxygen
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
concept bool semi_alphabet_concept = requires (t t1, t t2)
{
    // STL concepts
    requires std::is_pod_v<t> == true;
    requires std::is_swappable_v<t> == true;

    // static data members
    alphabet_size<t>::value;
    alphabet_size_v<t>;

    // conversion to rank
    { to_rank(t1) } -> underlying_rank_t<t>;

    // assignment from rank
    { assign_rank(t1,  0) } -> t &;
    { assign_rank(t{}, 0) } -> t &&;

    // required comparison operators
    { t1 == t2 } -> bool;
    { t1 != t2 } -> bool;
    { t1 <  t2 } -> bool;
    { t1 >  t2 } -> bool;
    { t1 <= t2 } -> bool;
    { t1 >= t2 } -> bool;
};
//!\endcond

/*!\interface seqan3::alphabet_concept <>
 * \extends seqan3::semi_alphabet_concept
 * \brief The generic alphabet concept that covers most data types used in ranges.
 * \ingroup alphabet
 *
 * \todo document
 *
 * \par Concepts and doxygen
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that satisfy this concept are shown as "implementing this interface".
 */
//!\cond
template <typename t>
concept bool alphabet_concept = requires (t t1, t t2)
{
    requires semi_alphabet_concept<t>;

    // conversion to char
    { to_char(t1) } -> underlying_char_t<t>;
    { std::cout << t1 };

    // assignment from char
    { assign_char(t1,  0) } -> t &;
    { assign_char(t{}, 0) } -> t &&;
};
//!\endcond

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
