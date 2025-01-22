// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::policy_affine_gap_with_trace_recursion.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/pairwise/detail/policy_affine_gap_recursion.hpp>

namespace seqan3::detail
{

/*!\brief Implements the alignment recursion function for the alignment algorithm using affine gap costs with trace
 *        information.
 * \ingroup alignment_pairwise
 * \copydetails seqan3::detail::policy_affine_gap_recursion
 */
template <typename alignment_configuration_t>
class policy_affine_gap_with_trace_recursion : protected policy_affine_gap_recursion<alignment_configuration_t>
{
protected:
    //!\brief The type of the base policy.
    using base_t = policy_affine_gap_recursion<alignment_configuration_t>;

    // Import base types.
    using typename base_t::affine_score_tuple_t;
    using typename base_t::score_type;
    using typename base_t::traits_type;

    //!\brief The trace type to use.
    using trace_type = typename traits_type::trace_type;
    //!\brief The internal tuple storing the trace directions of an affine cell.
    using affine_trace_tuple_t = std::tuple<trace_type, trace_type, trace_type>;
    //!\brief The affine cell type returned by the functions.
    using affine_cell_type = affine_cell_proxy<std::pair<affine_score_tuple_t, affine_trace_tuple_t>>;

    // Import base member.
    using base_t::first_column_is_free;
    using base_t::first_row_is_free;
    using base_t::gap_extension_score;
    using base_t::gap_open_score;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_affine_gap_with_trace_recursion() = default;                                               //!< Defaulted.
    policy_affine_gap_with_trace_recursion(policy_affine_gap_with_trace_recursion const &) = default; //!< Defaulted.
    policy_affine_gap_with_trace_recursion(policy_affine_gap_with_trace_recursion &&) = default;      //!< Defaulted.
    policy_affine_gap_with_trace_recursion &
    operator=(policy_affine_gap_with_trace_recursion const &) = default; //!< Defaulted.
    policy_affine_gap_with_trace_recursion &
    operator=(policy_affine_gap_with_trace_recursion &&) = default; //!< Defaulted.
    ~policy_affine_gap_with_trace_recursion() = default;            //!< Defaulted.

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::policy_affine_gap_recursion
    explicit policy_affine_gap_with_trace_recursion(alignment_configuration_t const & config) : base_t{config}
    {}
    //!\}

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::compute_inner_cell
    template <typename affine_cell_t>
    affine_cell_type compute_inner_cell(score_type diagonal_score,
                                        affine_cell_t previous_cell,
                                        score_type const sequence_score) const noexcept
    {
        diagonal_score += sequence_score;
        score_type horizontal_score = previous_cell.horizontal_score();
        score_type vertical_score = previous_cell.vertical_score();
        trace_directions best_trace = trace_directions::diagonal;

        diagonal_score = (diagonal_score < vertical_score)
                           ? (best_trace = previous_cell.vertical_trace(), vertical_score)
                           : (best_trace |= previous_cell.vertical_trace(), diagonal_score);
        diagonal_score =
            (diagonal_score < horizontal_score)
                ? (best_trace = previous_cell.horizontal_trace() | (best_trace & trace_directions::carry_up_open),
                   horizontal_score)
                : (best_trace |= previous_cell.horizontal_trace(), diagonal_score);

        score_type tmp = diagonal_score + gap_open_score;
        vertical_score += gap_extension_score;
        horizontal_score += gap_extension_score;

        // store the vertical_score and horizontal_score value in the next path
        trace_directions next_vertical_trace = trace_directions::up;
        trace_directions next_horizontal_trace = trace_directions::left;

        vertical_score =
            (vertical_score < tmp) ? (next_vertical_trace = trace_directions::up_open, tmp) : vertical_score;
        horizontal_score =
            (horizontal_score < tmp) ? (next_horizontal_trace = trace_directions::left_open, tmp) : horizontal_score;

        return {{diagonal_score, horizontal_score, vertical_score},
                {best_trace, next_horizontal_trace, next_vertical_trace}};
    }

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::initialise_origin_cell
    affine_cell_type initialise_origin_cell() const noexcept
    {
        return {base_t::initialise_origin_cell(),
                {trace_directions::none,
                 first_row_is_free ? trace_directions::none : trace_directions::left_open,
                 first_column_is_free ? trace_directions::none : trace_directions::up_open}};
    }

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::initialise_first_column_cell
    template <typename affine_cell_t>
    affine_cell_type initialise_first_column_cell(affine_cell_t previous_cell) const noexcept
    {
        return {base_t::initialise_first_column_cell(previous_cell),
                {previous_cell.vertical_trace(),
                 trace_directions::left_open,
                 first_column_is_free ? trace_directions::none : trace_directions::up}};
    }

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::initialise_first_row_cell
    template <typename affine_cell_t>
    affine_cell_type initialise_first_row_cell(affine_cell_t previous_cell) const noexcept
    {
        return {base_t::initialise_first_row_cell(previous_cell),
                {previous_cell.horizontal_trace(),
                 first_row_is_free ? trace_directions::none : trace_directions::left,
                 trace_directions::up_open}};
    }
};
} // namespace seqan3::detail
