// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::alignment_algorithm_state.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/matrix/detail/alignment_optimum.hpp>
#include <seqan3/utility/concept.hpp>
#include <seqan3/utility/simd/concept.hpp>

namespace seqan3::detail
{
/*!\brief Local state for the standard alignment algorithm.
 * \tparam score_type The type of the score; must model seqan3::arithmetic or seqan3::simd::simd_concept.
 * \ingroup alignment_pairwise
 *
 * \details
 *
 * This state is used internally for the standard alignment algorithm and caches the gap extension and gap open scores
 * as well as the current alignment optimum.
 * The alignment optimum stores the current score and the corresponding matrix coordinate in
 * the underlying two-dimensional matrix.
 */
template <typename score_type>
    requires arithmetic<score_type> || simd_concept<score_type>
struct alignment_algorithm_state
{
    //!\brief The cached gap extension score.
    score_type gap_extension_score{};
    //!\brief The cached gap open score.
    score_type gap_open_score{};
    //!\brief The current alignment optimum.
    alignment_optimum<score_type> optimum{};

    //!\brief Resets the alignment optimum to the default initialised optimum.
    constexpr void reset_optimum() noexcept
    {
        optimum = alignment_optimum<score_type>{};
    }
};

/*!\name Type deduction guides
 * \relates seqan3::detail::alignment_algorithm_state
 * \{
 */

//!\brief Deduces the template parameter for the score type from construction with gap open and gap extension scores.
template <typename score_type>
alignment_algorithm_state(score_type, score_type) -> alignment_algorithm_state<score_type>;
//!\}
} // namespace seqan3::detail
