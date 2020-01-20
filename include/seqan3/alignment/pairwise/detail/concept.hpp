// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concepts needed internally for the alignment algorithms.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\interface seqan3::detail::sequence_pair <>
 * \brief A helper concept to check if a type is a sequence pair.
 * \ingroup pairwise_alignment
 *
 * \tparam t The type to check.
 *
 * \details
 *
 * This concept checks if the given type models seqan3::tuple_like with exactly two elements and that both
 * types in the tuple model std::ranges::forward_range. Furthermore, the value type of both ranges must model
 * seqan3::semialphabet.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT sequence_pair = requires ()
{
    requires tuple_like<t>;
    requires std::tuple_size_v<t> == 2;
    requires std::ranges::forward_range<std::tuple_element_t<0, t>>;
    requires std::ranges::forward_range<std::tuple_element_t<1, t>>;
    requires semialphabet<std::ranges::range_value_t<std::tuple_element_t<0, t>>>;
    requires semialphabet<std::ranges::range_value_t<std::tuple_element_t<1, t>>>;
};
//!\endcond

/*!\interface seqan3::detail::sequence_pair_range <>
 * \brief A helper concept to check if a type is a range over seqan3::detail::sequence_pair.
 * \ingroup pairwise_alignment
 *
 * \tparam t The type to check.
 *
 * \details
 *
 * This concept checks if the given type models a std::ranges::forward_range and that the value type of the
 * range models seqan3::detail::sequence_pair.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT sequence_pair_range = std::ranges::forward_range<t> && sequence_pair<std::ranges::range_value_t<t>>;
//!\endcond

/*!\interface seqan3::detail::indexed_sequence_pair_range <>
 * \brief A helper concept to check the input of the range based alignment algorithm interface.
 * \ingroup pairwise_alignment
 *
 * \tparam t The type to check.
 *
 * \details
 *
 * This concept checks if the given type models a std::ranges::forward_range over indexed sequence pairs that are
 * passed to the alignment algorithms. An indexed sequence pair consists of a seqan3::detail::sequence_pair
 * that shall be aligned and an index that is used to identify the aligned sequence pair.
 * The caller can then infer the aligned sequences from the returned seqan3::alignment_result.
 * The layout of this indexed sequence type looks as follows:
 * * the first type of the pair must model seqan3::detail::sequence_pair, and
 * * the second type of the pair refers to the respective index type, which can be any type but must model
 *   std::copy_constructible.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT indexed_sequence_pair_range = std::ranges::forward_range<t> &&
                                             requires (std::ranges::range_value_t<t> value)
{
    requires tuple_like<decltype(value)>;
    requires std::tuple_size_v<decltype(value)> == 2;
    requires sequence_pair<std::tuple_element_t<0, decltype(value)>>;
    requires std::copy_constructible<std::tuple_element_t<1, decltype(value)>>;
};
//!\endcond

/*!\interface seqan3::detail::align_pairwise_single_input <>
 * \brief A helper concept to test for correct single value input in seqan3::align_pairwise.
 * \ingroup pairwise_alignment
 *
 * \tparam t The type to check.
 *
 * \details
 *
 * The given type must model seqan3::detail::sequence_pair and both contained types must model
 * std::ranges::viewable_range.
 *
 * \see seqan3::detail::align_pairwise_range_input
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT align_pairwise_single_input =
    sequence_pair<t> &&
    std::ranges::viewable_range<std::tuple_element_t<0, t>> &&
    std::ranges::viewable_range<std::tuple_element_t<1, t>>;
//!\endcond

/*!\interface seqan3::detail::align_pairwise_range_input <>
 * \brief A helper concept to test for correct range input in seqan3::align_pairwise.
 * \ingroup pairwise_alignment
 *
 * \tparam t The type to check.
 *
 * \details
 *
 * Only use input ranges whose value type models seqan3::detail::align_pairwise_value and
 * whose reference type is an lvalue reference and the range itself models std::ranges::viewable_range or
 * the reference type is a prvalue and it models seqan3::detail::align_pairwise_single_input.
 * This covers all typical use cases:
 * a) An lvalue range, whose reference type is a tuple like lvalue reference,
 * b) A range, whose reference type is a tuple over viewable ranges.
 * This covers also transforming and non-transforming views (e.g. views::zip, or views::take).
 * Only a temporary non-view range piped with views::persist can't be handled securely.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT align_pairwise_range_input =
    std::ranges::forward_range<t> &&
    sequence_pair<std::ranges::range_value_t<t>> &&
    ((std::ranges::viewable_range<t> && std::is_lvalue_reference_v<std::ranges::range_reference_t<t>>) ||
     align_pairwise_single_input<std::remove_reference_t<std::ranges::range_reference_t<t>>>);
//!\endcond
}  // namespace seqan3::detail
