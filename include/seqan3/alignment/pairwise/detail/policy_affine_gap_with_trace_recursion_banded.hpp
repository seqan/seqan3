// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_affine_gap_with_trace_recursion_banded.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/pairwise/detail/policy_affine_gap_with_trace_recursion.hpp>

namespace seqan3::detail
{

/*!\brief Implements the alignment recursion function for the banded alignment algorithm using affine gap costs with
 *        trace information.
 * \ingroup pairwise_alignment
 * \copydetails seqan3::detail::policy_affine_gap_recursion
 */
template <typename alignment_configuration_t>
class policy_affine_gap_with_trace_recursion_banded :
    protected policy_affine_gap_with_trace_recursion<alignment_configuration_t>
{
protected:
    //!\brief The type of the base policy.
    using base_t = policy_affine_gap_with_trace_recursion<alignment_configuration_t>;
    // Import base types.
    using typename base_t::traits_type;
    using typename base_t::score_type;
    using typename base_t::affine_cell_type;

    // Import base member.
    using base_t::gap_extension_score;
    using base_t::gap_open_score;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_affine_gap_with_trace_recursion_banded() = default; //!< Defaulted.
    //!\brief Defaulted.
    policy_affine_gap_with_trace_recursion_banded(policy_affine_gap_with_trace_recursion_banded const &) = default;
    //!\brief Defaulted.
    policy_affine_gap_with_trace_recursion_banded(policy_affine_gap_with_trace_recursion_banded &&) = default;
    policy_affine_gap_with_trace_recursion_banded & operator=(policy_affine_gap_with_trace_recursion_banded const &)
        = default; //!< Defaulted.
    policy_affine_gap_with_trace_recursion_banded & operator=(policy_affine_gap_with_trace_recursion_banded &&)
        = default; //!< Defaulted.
    ~policy_affine_gap_with_trace_recursion_banded() = default; //!< Defaulted.

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::policy_affine_gap_recursion
    explicit policy_affine_gap_with_trace_recursion_banded(alignment_configuration_t const & config) : base_t{config}
    {}
    //!\}

    //!\copydoc seqan3::detail::policy_affine_gap_recursion_banded::initialise_band_first_cell
    template <typename affine_cell_t>
    affine_cell_type initialise_band_first_cell(score_type diagonal_score,
                                                affine_cell_t previous_cell,
                                                score_type const sequence_score) const noexcept
    {
        diagonal_score += sequence_score;
        score_type horizontal_score = previous_cell.horizontal_score();
        trace_directions best_trace{};

        best_trace = previous_cell.horizontal_trace();
        diagonal_score = (diagonal_score < horizontal_score)
                       ? horizontal_score
                       : (best_trace |= trace_directions::diagonal, diagonal_score);

        score_type from_optimal_score = diagonal_score + gap_open_score;
        trace_directions next_horizontal_trace = trace_directions::left;

        horizontal_score += gap_extension_score;
        horizontal_score = (horizontal_score < from_optimal_score)
                         ? (next_horizontal_trace = trace_directions::left_open, from_optimal_score)
                         : horizontal_score;

        return {{diagonal_score, horizontal_score, from_optimal_score},
                {best_trace, next_horizontal_trace, trace_directions::up_open}};
    }
};
} // namespace seqan3::detail
