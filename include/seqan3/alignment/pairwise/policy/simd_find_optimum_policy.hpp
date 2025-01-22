// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::simd_find_optimum_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/alignment/matrix/detail/alignment_optimum.hpp>
#include <seqan3/alignment/pairwise/detail/alignment_algorithm_state.hpp>
#include <seqan3/alignment/pairwise/policy/find_optimum_policy.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/views/zip.hpp>

namespace seqan3::detail
{

/*!\brief A state that is only used for global alignments.
 * \ingroup alignment_pairwise
 * \tparam simd_t The simd vector score type; must model seqan3::simd::simd_concept.
 *
 * \details
 *
 * This state is only used for the global alignment to compute the correct optimum for every sequence pair.
 * If the sequences have different lengths the respective cells in the alignment matrix need to be queried to check for
 * the global optimum. In addition, the final scores and coordinates must be corrected as they are based on the
 * outer matrix defined by the longest sequence in the first and second collection.
 */
template <simd::simd_concept simd_t>
struct simd_global_alignment_state
{
    //!\brief The score offset that needs to be subtracted for every alignment to get the correct result.
    simd_t score_offset{};
    //!\brief A coordinate offset that needs to be subtracted for every alignment to get the correct end position.
    simd_t coordinate_offset{};
    //!\brief A mask vector storing the row indices for alignments that end in the last column of the global matrix.
    simd_t last_column_mask{};
    //!\brief A mask vector storing the column indices for alignments that end in the last row of the global matrix.
    simd_t last_row_mask{};
};

/*!\brief The CRTP-policy to determine the optimum of the dynamic programming matrix.
 * \ingroup alignment_pairwise_policy
 * \tparam alignment_algorithm_t The derived type (seqan3::detail::alignment_algorithm) to be augmented with this
 *                               CRTP-policy.
 * \tparam simd_t The simd vector type used for computing the scores; must model seqan3::simd::simd_concept.
 *
 * \details
 *
 * This class determines the matrix wide optimum. The search space can be further refined using the
 * `traits_type` which configures the search space of the alignment matrix.
 *
 * This class inherits from seqan3::detail::simd_global_alignment_state if it is a global alignment, otherwise it
 * inherits from seqan3::detail::empty_type that is optimised away due to empty base class optimisation.
 * The tracking behaviour for the optimum of the global alignment is slightly modified compared to the local alignment
 * or free end gap computation. Specifically, instead of checking for a new optimal score the coordinates are checked.
 * For the global alignment the cell which contains the optimum is already fixed (the sink of each matrix). To account
 * for different sizes of the sequences a padding match score is used to always add match scores after the end of a
 * smaller sequence has been reached. Accordingly, the target score will be mapped to a cell on
 * the last row or last column of the global matrix. Note the dimensions of the global matrix is defined by the size of
 * the longest sequence in collection 1 and 2 respectively. If the last cell of any of the contained matrices is
 * projected to a cell on the last row (following the diagonal going through the sink of this particular matrix),
 * then the column indices are tested to check if any sequence ends in this projected end point. Similarly, if the
 * last cell is projected to a cell on the last column, then the row indices are compared. The found score as well
 * as the indices are finally corrected to represent the original score and coordinates as if the sequence pair
 * was computed in scalar mode.
 */
template <typename alignment_algorithm_t, simd::simd_concept simd_t>
class simd_find_optimum_policy : public simd_global_alignment_state<simd_t>
{
private:
    //!\brief Befriends the derived class to grant it access to the private members.
    friend alignment_algorithm_t;

    //!\brief A flag to check if global alignment is computed.
    bool is_global_alignment{false};
    //!\brief Whether every cell of the alignment matrix shall be tracked.
    bool test_every_cell{false};
    //!\brief Whether cells of the last row shall be tracked.
    bool test_last_row_cell{false};
    //!\brief Whether cells of the last column shall be tracked.
    bool test_last_column_cell{false};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr simd_find_optimum_policy() = default;                                             //!< Defaulted.
    constexpr simd_find_optimum_policy(simd_find_optimum_policy const &) = default;             //!< Defaulted.
    constexpr simd_find_optimum_policy(simd_find_optimum_policy &&) = default;                  //!< Defaulted.
    constexpr simd_find_optimum_policy & operator=(simd_find_optimum_policy const &) = default; //!< Defaulted.
    constexpr simd_find_optimum_policy & operator=(simd_find_optimum_policy &&) = default;      //!< Defaulted.
    ~simd_find_optimum_policy() = default;                                                      //!< Defaulted.

    //!\brief Initialise the policy.
    template <typename configuration_t>
    simd_find_optimum_policy(configuration_t const & config)
    {
        if constexpr (configuration_t::template exists<align_cfg::method_local>())
            test_every_cell = true;

        is_global_alignment = configuration_t::template exists<align_cfg::method_global>();

        auto method_global_config = config.get_or(align_cfg::method_global{});

        test_last_row_cell = method_global_config.free_end_gaps_sequence1_trailing || is_global_alignment;
        test_last_column_cell = method_global_config.free_end_gaps_sequence2_trailing || is_global_alignment;
    }
    //!\}

protected:
    //!\copydoc seqan3::detail::find_optimum_policy::check_score_of_cell
    template <typename cell_t, typename score_t>
    constexpr void check_score_of_cell([[maybe_unused]] cell_t const & current_cell,
                                       [[maybe_unused]] alignment_algorithm_state<score_t> & state) const noexcept
    {
        if (test_every_cell)
            check_and_update(current_cell, state);
    }

    //!\brief Befriend the seqan3::detail::simd_affine_gap_policy to grant access to the check_score_of_cell function.
    template <typename other_alignment_algorithm_t, simd_concept score_t, typename is_local_t>
    friend class simd_affine_gap_policy;

    //!\brief Allow seqan3::detail::affine_gap_init_policy to access check_score.
    template <typename other_alignment_algorithm_t>
    friend class affine_gap_init_policy;

    //!\copydoc seqan3::detail::find_optimum_policy::check_score_of_last_row_cell
    template <typename cell_t, typename score_t>
    constexpr void
    check_score_of_last_row_cell([[maybe_unused]] cell_t const & last_row_cell,
                                 [[maybe_unused]] alignment_algorithm_state<score_t> & state) const noexcept
    {
        // Only search in last row if requested and not done already.}
        if (!test_every_cell && test_last_row_cell)
        {
            if (is_global_alignment)
                check_and_update<true>(last_row_cell, state);
            else
                check_and_update(last_row_cell, state);
        }
    }

    //!\copydoc seqan3::detail::find_optimum_policy::check_score_of_cells_in_last_column
    template <typename alignment_column_t, typename score_t>
    constexpr void
    check_score_of_cells_in_last_column([[maybe_unused]] alignment_column_t && last_column,
                                        [[maybe_unused]] alignment_algorithm_state<score_t> & state) const noexcept
    {
        // Only check last cell if not done before.
        if (!test_every_cell && test_last_column_cell)
        {
            for (auto && cell : last_column)
            {
                if (is_global_alignment)
                    check_and_update<false>(cell, state);
                else
                    check_and_update(cell, state);
            }
        }
    }

    //!\copydoc seqan3::detail::find_optimum_policy::check_score_of_last_cell
    template <typename cell_t, typename score_t>
    constexpr void check_score_of_last_cell([[maybe_unused]] cell_t const & last_cell,
                                            [[maybe_unused]] alignment_algorithm_state<score_t> & state) const noexcept
    {
        // Only check last cell if not done before.
        if (!(test_every_cell || test_last_row_cell || test_last_column_cell))
            check_and_update(last_cell, state);
    }

    /*!\brief Initialises the global alignment state for the current batch of sequences.
     *
     * \tparam sequence1_collection_t The type of the first collection; must model std::ranges::forward_range and
     *                                the value type must model std::ranges::forward_range.
     * \tparam sequence2_collection_t The type of the second collection; must model std::ranges::forward_range and
     *                                the value type must model std::ranges::forward_range.
     * \tparam score_t The type of the scoring scheme's padding score; must model seqan3::arithmetic.
     *
     * \param[in] sequence1_collection The first collection used for initialisation.
     * \param[in] sequence2_collection The second collection used for initialisation.
     * \param[in] padding_score The padding match score used for determining the resulting score offset.
     */
    template <std::ranges::forward_range sequence1_collection_t,
              std::ranges::forward_range sequence2_collection_t,
              arithmetic score_t>
    void initialise_find_optimum_policy([[maybe_unused]] sequence1_collection_t && sequence1_collection,
                                        [[maybe_unused]] sequence2_collection_t && sequence2_collection,
                                        [[maybe_unused]] score_t const padding_score)
    {
        if (is_global_alignment)
        {
            assert(std::ranges::distance(sequence1_collection) == std::ranges::distance(sequence2_collection));

            constexpr size_t simd_size = simd_traits<simd_t>::length;
            // First get global size.
            std::array<size_t, simd_size> sequence1_sizes{};
            std::array<size_t, simd_size> sequence2_sizes{};

            std::ptrdiff_t max_sequence1_size{};
            std::ptrdiff_t max_sequence2_size{};

            size_t array_index{};
            for (auto && [sequence1, sequence2] : views::zip(sequence1_collection, sequence2_collection))
            {
                sequence1_sizes[array_index] = std::ranges::distance(sequence1);
                sequence2_sizes[array_index] = std::ranges::distance(sequence2);
                max_sequence1_size = std::max<std::ptrdiff_t>(sequence1_sizes[array_index], max_sequence1_size);
                max_sequence2_size = std::max<std::ptrdiff_t>(sequence2_sizes[array_index], max_sequence2_size);
                ++array_index;
            }

            // The global diagonal ending in the sink of the outer alignment matrix.
            std::ptrdiff_t global_diagonal = max_sequence1_size - max_sequence2_size;

            for (size_t simd_index = 0; simd_index < simd_size; ++simd_index)
            {
                if (std::ptrdiff_t local_diagonal = sequence1_sizes[simd_index] - sequence2_sizes[simd_index];
                    local_diagonal < global_diagonal)
                { // optimum is stored in last row.
                    this->last_row_mask[simd_index] = max_sequence1_size - (global_diagonal - local_diagonal);
                    this->last_column_mask[simd_index] = max_sequence2_size + 1;
                    this->coordinate_offset[simd_index] = max_sequence2_size - sequence2_sizes[simd_index];
                    this->score_offset[simd_index] = padding_score * this->coordinate_offset[simd_index];
                }
                else // optimum is stored in last column
                {
                    this->last_column_mask[simd_index] = max_sequence2_size - (local_diagonal - global_diagonal);
                    this->last_row_mask[simd_index] = max_sequence1_size + 1;
                    this->coordinate_offset[simd_index] = max_sequence1_size - sequence1_sizes[simd_index];
                    this->score_offset[simd_index] = padding_score * this->coordinate_offset[simd_index];
                }
            }
        }
        // else no-op
    }

private:
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
        assert(!is_global_alignment); // This function should not be called for the global alignment.

        auto const & [score_cell, trace_cell] = cell;
        state.optimum.update_if_new_optimal_score(score_cell.current,
                                                  column_index_type{trace_cell.coordinate.first},
                                                  row_index_type{trace_cell.coordinate.second});
    }

    /*!\brief Tests if the current row, respectively column, is part of a global alignment to track.
     *
     * \tparam in_last_row A bool constant to indicate whether this check is done for the last row.
     * \tparam cell_t The type of the alignment matrix cell. The cell type corresponds to the value type of the range
     *                returned by seqan3::detail::alignment_matrix_policy::current_alignment_column.
     * \tparam score_t The alignment algorithm score type.
     *
     * \param[in] cell The current cell to get the score and the coordinate from.
     * \param[in,out] state The state with the current optimum to update.
     *
     * \details
     *
     * Checks for the last row if the current column coordinate matches any column coordinate stored in the mask for
     * the last row, respectively for the last column if the current row coordinate matches any row coordinate stored
     * in the mask for the last column.
     */
    template <bool in_last_row, typename cell_t, typename score_t>
    constexpr void check_and_update(cell_t const & cell, alignment_algorithm_state<score_t> & state) const noexcept
    {
        assert(is_global_alignment); // This function should only be called for the global alignment.

        using simd_mask_t = typename simd_traits<simd_t>::mask_type;
        auto const & [score_cell, trace_cell] = cell;
        simd_t column_positions = simd::fill<simd_t>(trace_cell.coordinate.first);
        simd_t row_positions = simd::fill<simd_t>(trace_cell.coordinate.second);

        simd_mask_t mask{};

        if constexpr (in_last_row) // Check if column was masked
            mask = (column_positions == this->last_row_mask);
        else // Check if row was masked
            mask = (row_positions == this->last_column_mask);

        // In global alignment we are only interested in the position not the max of the scores.
        // In addition, the scores need to be corrected in order to track the right score.
        state.optimum.score = mask ? score_cell.current - this->score_offset : state.optimum.score;
        state.optimum.column_index = mask ? column_positions - this->coordinate_offset : state.optimum.column_index;
        state.optimum.row_index = mask ? row_positions - this->coordinate_offset : state.optimum.row_index;
    }
};

} // namespace seqan3::detail
