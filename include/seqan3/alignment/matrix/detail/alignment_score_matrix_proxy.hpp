// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_score_matrix_proxy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/simd/concept.hpp>

namespace seqan3::detail
{

/*!\brief A proxy type for a unified access to the score matrix during alignment computation.
 * \tparam score_type The type wrapped by this proxy.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * Provides named accessors to the respective values of the score matrix during the alignment computation.
 * The `score_type` must either model the seqan3::arithmetic or seqan3::detail::simd_conceptvector concept.
 * In some cases the last value to the left is stored in a different location then the next left value that
 * will be used in the subsequent column (e.g. in the banded alignment matrix). Accordingly the left value is
 * split in a reference for reading, referred to as `r_left` (for reading), and a reference to the next left value
 * on the same row, referred to as `w_left` (for writing).
 */
template <typename score_type>
struct alignment_score_matrix_proxy
{
    static_assert(arithmetic<score_type> || simd_concept<score_type>,
                  "Value type must either be an arithmetic type or a simd vector type.");

    score_type & current; //!< The current value.
    score_type & diagonal; //!< The last diagonal value.
    score_type & r_left; //!< The last value to the left.
    score_type & w_left; //!< The next value to the left.
    score_type & up; //!< The last value above.
};
} // namespace seqan3::detail
