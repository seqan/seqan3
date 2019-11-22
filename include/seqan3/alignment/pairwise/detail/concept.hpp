// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
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
 * This concept checks if the given type models a seqan3::tuple_like type with exactly two elements and that both
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
 * \brief A helper concept to check if a type is a range over seqan3::detail::sequence_pair's.
 * \ingroup pairwise_alignment
 *
 * \tparam t The type to check.
 *
 * \details
 *
 * This concept checks if the given type models a std::ranges::forward_range over and that the value type of the
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
 * * the second type of the pair refers to the respective index type.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT indexed_sequence_pair_range = std::ranges::forward_range<t> &&
                                             requires (std::ranges::range_value_t<t> value)
{
    requires tuple_like<decltype(value)>;
    requires std::tuple_size_v<decltype(value)> == 2;
    requires sequence_pair<std::tuple_element_t<0, decltype(value)>>;
};
//!\endcond
}  // namespace seqan3::detail
