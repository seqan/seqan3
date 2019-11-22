// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides utility traits for the configuration and execution of the alignment algorithm.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/pairwise/detail/concept.hpp>
#include <seqan3/range/views/chunk.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief A transformation trait to retrieve the chunked range over indexed sequence pairs.
 * \ingroup pairwise_alignment
 * \implements seqan3::transformation_trait
 *
 * \tparam sequence_pairs_t The type of the sequences to be transformed; must model seqan3::detail::sequence_pair_range.
 *
 * \details
 *
 * This transformation trait transforms a range over sequence pairs into an indexed range over sequence pairs.
 * In addition, the range is chunked which is the common interface for alignment algorithms.
 * The returned type models seqan3::detail::indexed_sequence_pair_range.
 */
template <typename sequence_pairs_t>
//!\cond
    requires sequence_pair_range<std::remove_reference_t<sequence_pairs_t>>
//!\endcond
struct chunked_indexed_sequence_pairs
{
    //!\brief The transformed type that models seqan3::detail::indexed_sequence_pair_range.
    using type = decltype(views::zip(std::declval<sequence_pairs_t>(), std::views::iota(0)) | views::chunk(1));
};

}  // namespace seqan3::detail
