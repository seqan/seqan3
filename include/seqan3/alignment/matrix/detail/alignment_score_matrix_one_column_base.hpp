// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_score_matrix_one_column_base.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <utility>
#include <vector>

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>

namespace seqan3::detail
{

/*!\brief A base class for alignment score matrices using only one column to compute the matrix.
 * \tparam score_t The type of the scores.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * Manages the actual storage as a std::vector. How much memory is allocated is handled by the derived type.
 * The `score_t` must either model the seqan3::arithmetic or seqan3::detail::simd_conceptvector concept.
 */
template <typename score_t>
struct alignment_score_matrix_one_column_base
{
protected:
    static_assert(arithmetic<score_t> || simd_concept<score_t>,
                  "Score type must either be either an arithmetic type or a simd vector type.");

    //!\brief The underlying type of the scores.
    using underlying_type = score_t;
    //!\brief The actual element type.
    using element_type = std::tuple<underlying_type, underlying_type>;
    //!\brief The allocator type.
    using allocator_type = aligned_allocator<element_type, sizeof(element_type)>;
    //!\brief The type of the underlying storage.
    using pool_type = std::vector<element_type, allocator_type>;
    //!\brief The size type.
    using size_type = size_t;
public:
    //!\brief The linearised memory pool storing only one column of the matrix.
    pool_type pool{};
    //!\brief Internal cache for the last diagonal and vertical value during the alignment computation.
    std::array<underlying_type, 3> cache{};  // Third argument is used to cache next diagonal value in non-banded case.
    //!\brief The number of columns.
    size_type num_cols{};
    //!\brief The number of num_rows.
    size_type num_rows{};
};

} // namespace seqan3::detail
