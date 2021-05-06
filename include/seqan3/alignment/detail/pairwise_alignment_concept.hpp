// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::pairwise_alignment and seqan3::detail::writable_pairwise_alignment.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#if SEQAN3_WORKAROUND_ISSUE_286
#include <seqan3/core/detail/transfer_type_modifier_onto.hpp>
#endif // SEQAN3_WORKAROUND_ISSUE_286
#include <seqan3/utility/tuple/concept.hpp>

namespace seqan3::detail
{

/*!\interface seqan3::detail::pairwise_alignment < >
 * \extends seqan3::pair_like
 * \ingroup alignment
 * \brief A concept that models a pairwise alignment type.
 *
 * \details
 *
 * A pairwise alignment is a seqan3::pair_like of two seqan3::aligned_sequence's.
 */
//!\cond
template <typename pairwise_alignment_t>
SEQAN3_CONCEPT pairwise_alignment =
    pair_like<pairwise_alignment_t> &&
#if SEQAN3_WORKAROUND_ISSUE_286
    aligned_sequence<
        // simulate that tuple_element transfers const over to it's inner members.
        transfer_type_modifier_onto_t<
            std::remove_reference_t<pairwise_alignment_t>,
            std::tuple_element_t<0, std::remove_cvref_t<pairwise_alignment_t>>
        >
    > &&
    aligned_sequence<
        // simulate that tuple_element transfers const over to it's inner members.
        transfer_type_modifier_onto_t<
            std::remove_reference_t<pairwise_alignment_t>,
            std::tuple_element_t<1, std::remove_cvref_t<pairwise_alignment_t>>
        >
    >;
#else // ^^^ workaround / no workaround vvv
    aligned_sequence<std::tuple_element_t<0, std::remove_reference_t<pairwise_alignment_t>>> &&
    aligned_sequence<std::tuple_element_t<1, std::remove_reference_t<pairwise_alignment_t>>>;
#endif // SEQAN3_WORKAROUND_ISSUE_286

//!\endcond

/*!\interface seqan3::detail::writable_pairwise_alignment < >
 * \extends seqan3::detail::pairwise_alignment
 * \ingroup alignment
 * \brief A concept that models a writable pairwise alignment type.
 *
 * \details
 *
 * A writable pairwise alignment is a seqan3::pair_like of two seqan3::writable_aligned_sequence's.
 */
//!\cond
template <typename pairwise_alignment_t>
SEQAN3_CONCEPT writable_pairwise_alignment =
    pairwise_alignment<pairwise_alignment_t> &&
#if SEQAN3_WORKAROUND_ISSUE_286
    writable_aligned_sequence<
        // simulate that tuple_element transfers const over to it's inner members.
        transfer_type_modifier_onto_t<
            std::remove_reference_t<pairwise_alignment_t>,
            std::tuple_element_t<0, std::remove_cvref_t<pairwise_alignment_t>>
        >
    > &&
    writable_aligned_sequence<
        // simulate that tuple_element transfers const over to it's inner members.
        transfer_type_modifier_onto_t<
            std::remove_reference_t<pairwise_alignment_t>,
            std::tuple_element_t<1, std::remove_cvref_t<pairwise_alignment_t>>
        >
    >;
#else // ^^^ workaround / no workaround vvv
    writable_aligned_sequence<std::tuple_element_t<0, std::remove_reference_t<pairwise_alignment_t>>> &&
    writable_aligned_sequence<std::tuple_element_t<1, std::remove_reference_t<pairwise_alignment_t>>>;
#endif // SEQAN3_WORKAROUND_ISSUE_286
//!\endcond

} // namespace seqan3::detail
