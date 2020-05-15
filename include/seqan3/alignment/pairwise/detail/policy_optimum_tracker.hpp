// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_optimum_tracker.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <limits>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>

namespace seqan3::detail
{

/*!\brief Implements the tracker to store the global optimum for a particular alignment computation.
 * \ingroup pairwise_alignment
 *
 * \tparam alignment_configuration_t The type of the alignment configuration; must be a type specialisation of
 *                                   seqan3::configuration.
 *
 * \details
 *
 * Implements the interface to track the alignment optimum. It updates the currently stored optimum if the
 * optimal score of the given cell compares greater or equal to the stored optimum.
 * Special methods are offered to track any cell (for example when computing the local alignment), the last cell of a
 * column or a row (for example when using free-end gaps), or the final cell of the entire matrix (for example in the
 * standard global alignment).
 * Depending on the configuration they might or might not evaluate the stored optimum score of that cell to compare
 * it with the currently stored optimum.
 * The optimum needs to be reset in between alignment computations in order to ensure that the correct result is
 * tracked.
 */
template <typename alignment_configuration_t>
//!\cond
    requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
//!\endcond
class policy_optimum_tracker
{
private:
    //!\brief The configuration traits type.
    using traits_type = alignment_configuration_traits<alignment_configuration_t>;
    //!\brief The configured score type.
    using score_type = typename traits_type::score_type;

    //!\brief The tracked score of the global optimum.
    score_type optimal_score{};
    //!\brief The matrix coordinate of the tracked optimum.
    matrix_coordinate optimal_coordinate{};

protected:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_optimum_tracker() = default; //!< Defaulted.
    policy_optimum_tracker(policy_optimum_tracker const &) = default; //!< Defaulted.
    policy_optimum_tracker(policy_optimum_tracker &&) = default; //!< Defaulted.
    policy_optimum_tracker & operator=(policy_optimum_tracker const &) = default; //!< Defaulted.
    policy_optimum_tracker & operator=(policy_optimum_tracker &&) = default; //!< Defaulted.
    ~policy_optimum_tracker() = default; //!< Defaulted.

    /*!\brief Construction and initialisation using the alignment configuration.
     * \param[in] config The alignment configuration (not used in this context).
     *
     * \details
     *
     * Resets the optimum on construction.
     */
    policy_optimum_tracker(alignment_configuration_t const & SEQAN3_DOXYGEN_ONLY(config))
    {
        reset_optimum();
    }
    //!\}

    /*!\brief Tracks any cell within the alignment matrix.
     *
     * \tparam cell_t The cell type of the alignment matrix; must have a member function `optimal_score()`.
     *
     * \param[in] cell The current cell to be tracked.
     * \param[in] coordinate The matrix coordinate of the current cell.
     *
     * \returns The forwarded cell.
     *
     * \details
     *
     * A call to this function only tracks the optimal score of the given cell if the configuration of the alignment
     * algorithm requires it, for example when a local alignment shall be computed.
     */
    template <typename cell_t>
    decltype(auto) track_cell(cell_t && cell, matrix_coordinate SEQAN3_DOXYGEN_ONLY(coordinate)) noexcept
    {
        return std::forward<cell_t>(cell);
    }

    /*!\brief Tracks the last cell of a row within the alignment matrix.
     *
     * \tparam cell_t The cell type of the alignment matrix; must have a member function `optimal_score()`.
     *
     * \param[in] cell The current cell to be tracked.
     * \param[in] coordinate The matrix coordinate of the current cell.
     *
     * \returns The forwarded cell.
     *
     * \details
     *
     * A call to this function only tracks the optimal score of the given cell if the configuration of the alignment
     * algorithm requires it, for example when a semi-global alignment shall be computed.
     */
    template <typename cell_t>
    decltype(auto) track_last_row_cell(cell_t && cell, matrix_coordinate SEQAN3_DOXYGEN_ONLY(coordinate)) noexcept
    {
        return std::forward<cell_t>(cell);
    }

    /*!\brief Tracks the last cell of a column within the alignment matrix.
     *
     * \tparam cell_t The cell type of the alignment matrix; must have a member function `optimal_score()`.
     *
     * \param[in] cell The current cell to be tracked.
     * \param[in] coordinate The matrix coordinate of the current cell.
     *
     * \returns The forwarded cell.
     *
     * \details
     *
     * A call to this function only tracks the optimal score of the given cell if the configuration of the alignment
     * algorithm requires it, for example when a semi-global alignment shall be computed.
     */
    template <typename cell_t>
    decltype(auto) track_last_column_cell(cell_t && cell, matrix_coordinate SEQAN3_DOXYGEN_ONLY(coordinate)) noexcept
    {
        return std::forward<cell_t>(cell);
    }

    /*!\brief Tracks the final cell of the alignment matrix.
     *
     * \tparam cell_t The cell type of the alignment matrix; must have a member function `optimal_score()`.
     *
     * \param[in] cell The current cell to be tracked.
     * \param[in] coordinate The matrix coordinate of the current cell.
     *
     * \returns The forwarded cell.
     *
     * \details
     *
     * A call to this function only tracks the optimal score of the given cell if the configuration of the alignment
     * algorithm requires it, for example when a global alignment shall be computed.
     */
    template <typename cell_t>
    decltype(auto) track_final_cell(cell_t && cell, matrix_coordinate coordinate) noexcept
    {
        bool const is_better_score = cell.optimal_score() >= optimal_score;
        optimal_score = (is_better_score) ? cell.optimal_score() : optimal_score;
        optimal_coordinate = (is_better_score) ? std::move(coordinate) : optimal_coordinate;

        return std::forward<cell_t>(cell);
    }

    /*!\brief Returns the stored optimum.
     *
     * \returns The stored optimum.
     */
    score_type tracked_optimum() const noexcept
    {
        return optimal_score;
    }

    //!\brief Resets the optimum such that a new alignment can be computed.
    void reset_optimum() noexcept
    {
        optimal_score = std::numeric_limits<score_type>::lowest();
        optimal_coordinate = {};
    }
};
} // namespace seqan3::detail
