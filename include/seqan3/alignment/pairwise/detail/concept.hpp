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

#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief A helper concept to check the input of the range based alignment algorithm interface.
 * \ingroup pairwise_alignment
 *
 * \tparam t The type to check.
 *
 * \details
 *
 * This concept checks if the given type models a std::ranges::forward_range over indexed sequence pairs that are
 * passed to the alignment algorithms. An indexed sequence pair is a pair of sequences that shall be aligned. In
 * addition, the pair is annotated with an index so that the caller can infer the aligned sequences from the returned
 * seqan3::alignment_result. The layout of this indexed sequence type looks as follows:
 * * the first type of the pair refers to the sequence pair type, and
 * * the second type of the pair refers to the respective index type.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT indexed_sequence_pair_range = requires (std::ranges::range_value_t<t> value)
{
    requires std::ranges::forward_range<t>;
    requires tuple_like<decltype(value)>;
    requires std::tuple_size_v<decltype(value)> == 2;
    requires tuple_like<std::tuple_element_t<0, decltype(value)>>;
    requires std::tuple_size_v<std::tuple_element_t<0, decltype(value)>> == 2;
};
//!\endcond
}  // namespace seqan3::detail
