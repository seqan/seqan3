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
 * \brief Includes the aligned_sequence_concept and the related insert_gap and
 *        erase_gap functions to enable stl container support.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/exception.hpp>
#include <seqan3/alphabet/gap/all.hpp>
#include <seqan3/core/metafunction/all.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{

// -----------------------------------------------------------------------------
// aligned_sequence_concept
// -----------------------------------------------------------------------------

/*!\interface seqan3::aligned_sequence_concept <>
 * \extends   std::ranges::ForwardRange
 * \brief     The generic concept for an aligned sequence.
 * \ingroup   aligned_sequence
 *
 * This concept describes the requirements a sequence must fulfil
 * in order to be part of the seqan3::alignment object.
 *
 * The following extended type requirements for a type `T` must hold true:
 *
 *   * seqan3::reference_t<T> must model seqan3::alphabet_concept.
 *   * seqan3::reference_t<T> must be assignable from seqan3::gap.
 *
 * ### Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and metafunctions.
 * Types that model this concept are shown as "implementing this interface".
 */
/*!\name Requirements for seqan3::aligned_sequence_concept
 * \brief You can expect these functions on all types that model seqan3::aligned_sequence_concept.
 * \relates seqan3::aligned_sequence_concept
 * \{
 */
/*!\fn      inline seq_type::iterator insert_gap(seq_type & seq, typename seq_type::const_iterator pos_it)
 * \brief   Insert a seqan3::gap into an aligned sequence.
 * \relates seqan3::aligned_sequence_concept
 *
 * \tparam        seq_type Type of the range to modify; must model
 *                         seqan3::aligned_sequence_concept.
 * \param[in,out] seq      The aligned sequence to modify.
 * \param[in]     pos_it   The iterator pointing to the position where to insert a gap.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn      inline seq_type::iterator insert_gap(seq_type & seq, typename seq_type::const_iterator pos_it, typename seq_type::size_type count)
 * \brief   Insert multiple seqan3::gap into an aligned sequence.
 * \relates seqan3::aligned_sequence_concept
 *
 * \tparam        seq_type Type of the range to modify; must model
 *                         seqan3::aligned_sequence_concept.
 * \param[in,out] seq      The aligned sequence to modify.
 * \param[in]     pos_it   The iterator pointing to the position where to insert a gaps.
 * \param[in]     count    The number of gaps to insert.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn      inline seq_type::iterator erase_gap(seq_type & seq, typename seq_type::const_iterator pos_it)
 * \brief   Erase a seqan3::gap from an aligned sequence.
 * \relates seqan3::aligned_sequence_concept
 *
 * \tparam        seq_type Type of the range to modify; must model
 *                         seqan3::aligned_sequence_concept.
 * \param[in,out] seq      The aligned sequence to modify.
 * \param[in]     pos_it   The iterator pointing to the position where to erase a gap.
 *
 * \throws seqan3::gap_erase_failure if there is no seqan3::gap at \p pos_it.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn      inline seq_type::iterator erase_gap(seq_type & seq, typename seq_type::const_iterator first, typename seq_type::const_iterator last)
 * \brief   Erase multiple seqan3::gap from an aligned sequence.
 * \relates seqan3::aligned_sequence_concept
 *
 * \tparam        seq_type Type of the range to modify; must model
 *                         seqan3::aligned_sequence_concept.
 * \param[in,out] seq      The aligned sequence to modify.
 * \param[in]     first    The iterator pointing to the position where to start erasing gaps.
 * \param[in]     last     The iterator pointing to the position where to stop erasing gaps.
 *
 * \throws seqan3::gap_erase_failure if one of the characters in [\p first, \p last) no seqan3::gap.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
//!\}
//!\cond
template <typename t>
concept aligned_sequence_concept =
    std::ranges::ForwardRange<t> &&
    alphabet_concept<value_type_t<t>> &&
    std::Assignable<reference_t<t>, gap const &> &&
    requires (t v)
    {
        { insert_gap(v, v.begin()) } -> typename t::iterator; // global functions for generic usability
        { insert_gap(v, v.begin(), 2) } -> typename t::iterator;
        { erase_gap(v, v.begin()) } -> typename t::iterator;
        { erase_gap(v, v.begin(), v.end()) } -> typename t::iterator;
    };
//!\endcond

} // namespace seqan
