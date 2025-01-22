// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::find_optimum_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/matrix/detail/alignment_optimum.hpp>
#include <seqan3/alignment/pairwise/detail/alignment_algorithm_state.hpp>

namespace seqan3::detail
{

/*!\brief The CRTP-policy to determine the optimum of the dynamic programming matrix.
 * \ingroup alignment_pairwise_policy
 * \tparam alignment_algorithm_t The derived type (seqan3::detail::alignment_algorithm) to be augmented with this
 *                               CRTP-policy.
 *
 * \details
 *
 * This class determines the matrix wide optimum. The search space can be further refined using the
 * `traits_type` which configures the search space of the alignment matrix.
 */
template <typename alignment_algorithm_t>
class find_optimum_policy
{
private:
    //!\brief Befriends the derived class to grant it access to the private members.
    friend alignment_algorithm_t;

    //!\brief Whether every cell of the alignment matrix shall be tracked.
    bool test_every_cell{false};
    //!\brief Whether cells of the last row shall be tracked.
    bool test_last_row_cell{false};
    //!\brief Whether cells of the last column shall be tracked.
    bool test_last_column_cell{false};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr find_optimum_policy() = default;                                        //!< Defaulted
    constexpr find_optimum_policy(find_optimum_policy const &) = default;             //!< Defaulted
    constexpr find_optimum_policy(find_optimum_policy &&) = default;                  //!< Defaulted
    constexpr find_optimum_policy & operator=(find_optimum_policy const &) = default; //!< Defaulted
    constexpr find_optimum_policy & operator=(find_optimum_policy &&) = default;      //!< Defaulted
    ~find_optimum_policy() = default;                                                 //!< Defaulted

    //!\brief Initialises the policy with the configuration.
    template <typename configuration_t>
    find_optimum_policy(configuration_t const & config)
    {
        if constexpr (configuration_t::template exists<align_cfg::method_local>())
            test_every_cell = true;

        auto method_global_config = config.get_or(align_cfg::method_global{});
        test_last_row_cell = method_global_config.free_end_gaps_sequence1_trailing;
        test_last_column_cell = method_global_config.free_end_gaps_sequence2_trailing;
    }
    //!\}

protected:
    /*!\brief Checks if a given cell is a new optimum in the alignment.
     * \tparam cell_t The type of the alignment matrix cell.
     * \tparam score_t The type of the score.
     *
     * \param[in] current_cell The currently computed alignment matrix cell.
     * \param[in,out] state The state with the current optimum to update.
     *
     * \details
     *
     * This function resolves to a "NO-OP" function if the trait for searching every cell is set to std::false_type.
     */
    template <typename cell_t, typename score_t>
    constexpr void check_score_of_cell([[maybe_unused]] cell_t const & current_cell,
                                       [[maybe_unused]] alignment_algorithm_state<score_t> & state) const noexcept
    {
        if (test_every_cell)
            check_and_update(current_cell, state);
    }

    //!\brief Allow seqan3::detail::affine_gap_policy to access check_score.
    template <typename other_alignment_algorithm_t, typename score_t, typename is_local_t>
    friend class affine_gap_policy;

    template <typename other_alignment_algorithm_t, simd_concept score_t, typename is_local_t>
    friend class simd_affine_gap_policy;

    //!\brief Allow seqan3::detail::affine_gap_init_policy to access check_score.
    template <typename other_alignment_algorithm_t>
    friend class affine_gap_init_policy;

    /*!\brief Checks if a cell in the last row of the alignment matrix is a new optimum in the alignment.
     * \tparam cell_t The type of the alignment matrix cell.
     * \tparam score_t The type of the score.
     *
     * \param[in] last_row_cell The cell of the current column in the last row of the alignment matrix.
     * \param[in,out] state The state with the current optimum to update.
     *
     * \details
     *
     * This function resolves to a "NO-OP" function if the trait for searching the last row is set to std::false_type
     * or if the trait for searching every cell is set to std::true_type.
     */
    template <typename cell_t, typename score_t>
    constexpr void
    check_score_of_last_row_cell([[maybe_unused]] cell_t const & last_row_cell,
                                 [[maybe_unused]] alignment_algorithm_state<score_t> & state) const noexcept
    {
        // Only search in last row if requested and not done already.
        if (!test_every_cell && test_last_row_cell)
            check_and_update(last_row_cell, state);
    }

    /*!\brief Checks all cells of the last alignment column for a new alignment optimum.
     * \tparam alignment_column_t The type of an alignment column.
     * \tparam score_t The type of the optimal score.
     *
     * \param[in] last_column The last column of the alignment matrix.
     * \param[in,out] state The state with the current optimum to update.
     *
     * \details
     *
     * This function resolves to a "NO-OP" function if the trait for searching the last column is set to std::false_type
     * or if the trait for searching every cell is set to std::true_type.
     */
    template <typename alignment_column_t, typename score_t>
    constexpr void
    check_score_of_cells_in_last_column([[maybe_unused]] alignment_column_t && last_column,
                                        [[maybe_unused]] alignment_algorithm_state<score_t> & state) const noexcept
    {
        // Only check last cell if not done before.
        if (!test_every_cell && test_last_column_cell)
            for (auto && cell : last_column)
                check_and_update(cell, state);
    }

    /*!\brief Checks if the last cell of the alignment matrix is a new optimum in the alignment.
     * \tparam cell_t The type of the last cell.
     * \tparam score_t The type of the score.
     *
     * \param[in] last_cell The last cell of the alignment matrix.
     * \param[in,out] state The state with the current optimum to update.
     *
     * \details
     *
     * This function resolves to a "NO-OP" function if the last cell has been checked already as part of the last
     * row, last column, or in case every cell was checked.
     */
    template <typename cell_t, typename score_t>
    constexpr void check_score_of_last_cell([[maybe_unused]] cell_t const & last_cell,
                                            [[maybe_unused]] alignment_algorithm_state<score_t> & state) const noexcept
    {
        // Only check last cell if not done before.
        if (!(test_every_cell || test_last_row_cell || test_last_column_cell))
            check_and_update(last_cell, state);
    }

    /*!\brief Tests if the score in the current cell is greater than the current alignment optimum.
     * \tparam cell_t The type of the alignment matrix cell. The cell type corresponds to the value type of the range
     *                returned by seqan3::detail::alignment_matrix_policy::current_alignment_column.
     * \tparam score_t The alignment algorithm score type.
     *
     * \param[in] cell The current cell to get the score and the coordinate from.
     * \param[in,out] state The state with the current optimum to update.
     */
    template <typename cell_t, typename score_t>
    constexpr void check_and_update(cell_t const & cell, alignment_algorithm_state<score_t> & state) const noexcept
    {
        auto const & [score_cell, trace_cell] = cell;
        state.optimum.update_if_new_optimal_score(score_cell.current,
                                                  column_index_type{trace_cell.coordinate.first},
                                                  row_index_type{trace_cell.coordinate.second});
    }
};

} // namespace seqan3::detail
