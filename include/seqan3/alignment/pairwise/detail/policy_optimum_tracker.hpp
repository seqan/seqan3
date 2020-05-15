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

#include <seqan3/alignment/matrix/detail/coordinate_matrix.hpp>
#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>

namespace seqan3::detail
{

/*!\brief A binary operation to update the alignment optimum.
 * \ingroup pairwise_alignment
 *
 * \details
 *
 * Updates the current alignment optimum with the new score and the respective coordinate if the new score
 * compares greater or equal to the score of the current optimum.
 */
struct alignment_optimum_updater_greater_equal
{
    /*!\brief The binary compare and update operation.
     * \tparam lhs_t The type of the left hand side; must model seqan3::tuple_like with a tuple size of 2.
     * \tparam rhs_t The type of the right hand side; must model seqan3::tuple_like with a tuple size of 2.
     *
     * \param[in,out] optimal_score_coordinate_pair The current optimum.
     * \param[in] current_cell_score_coordinate_pair The new score and coordinate to compare with.
     *
     * \details
     *
     * Requires that first value of the tuple represents the score and the second type the matrix coordinate.
     */
    template <tuple_like lhs_t, tuple_like rhs_t>
    //!\cond
        requires std::tuple_size_v<lhs_t> == 2 && std::tuple_size_v<rhs_t> == 2 &&
                 std::totally_ordered_with<std::tuple_element_t<0, lhs_t>, std::tuple_element_t<0, rhs_t>> &&
                 std::assignable_from<std::tuple_element_t<0, lhs_t>, std::tuple_element_t<0, rhs_t>> &&
                 std::assignable_from<std::tuple_element_t<1, lhs_t>, std::tuple_element_t<1, rhs_t>>
    //!\endcond
    void operator()(lhs_t && optimal_score_coordinate_pair,
                    rhs_t && current_cell_score_coordinate_pair) const
    {
        auto && [optimal_score, optimal_coordinate] = optimal_score_coordinate_pair;
        auto && [current_cell_score, current_cell_coordinate] = current_cell_score_coordinate_pair;

        bool const is_better_score = current_cell_score >= optimal_score;
        optimal_score = (is_better_score) ? current_cell_score : optimal_score;
        optimal_coordinate = (is_better_score) ? std::move(current_cell_coordinate) : optimal_coordinate;
    }
};

/*!\brief Implements the tracker to store the global optimum for a particular alignment computation.
 * \ingroup pairwise_alignment
 *
 * \tparam alignment_configuration_t The type of the alignment configuration; must be a type specialisation of
 *                                   seqan3::configuration.
 * \tparam binary_update_operation_t The type of a binary update operation to update the alignment optimum; must model
 *                                   std::semiregular and must be invocable with two tuples containing the score and
 *                                   the coordinate.
 * \details
 *
 * Implements the interface to track the alignment optimum. It updates the currently stored optimum using the
 * binary update operation. The binary update operation is stored inside of the class and can have a state.
 * Special methods are offered to track any cell (for example when computing the local alignment), the last cell of a
 * column or a row (for example when using free-end gaps), or the final cell of the entire matrix (for example in the
 * standard global alignment).
 * The optimum needs to be reset in between alignment computations in order to ensure that the correct result is
 * tracked.
 */
template <typename alignment_configuration_t, std::semiregular binary_update_operation_t>
//!\cond
    requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
//!\endcond
class policy_optimum_tracker
{
protected:
    //!\brief The configuration traits type.
    using traits_type = alignment_configuration_traits<alignment_configuration_t>;
    //!\brief The configured score type.
    using score_type = typename traits_type::score_type;
    //!\brief The matrix coordinate type that is used to locate a cell inside of the alignment matrix.
    using matrix_coordinate_type = typename traits_type::matrix_coordinate_type;

    /*!\name Data parameters
     * \{
     */
    //!\brief The tracked score of the global optimum.
    score_type optimal_score{};
    //!\brief The matrix coordinate of the tracked optimum.
    matrix_coordinate_type optimal_coordinate{};
    //!\brief The function object to compare and exchange the optimum.
    binary_update_operation_t binary_update_operation{};
    //!\}

    /*!\name Switch parameters
     * \{
     */
    //!\brief Whether every cell of the alignment matrix shall be tracked.
    bool test_every_cell{false};
    //!\brief Whether cells of the last row shall be tracked.
    bool test_last_row_cell{false};
    //!\brief Whether cells of the last column shall be tracked.
    bool test_last_column_cell{false};
    //!\}

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
    decltype(auto) track_cell(cell_t && cell, matrix_coordinate_type coordinate) noexcept
    {
        if (test_every_cell)
            binary_update_operation(std::forward_as_tuple(optimal_score, optimal_coordinate),
                                    std::forward_as_tuple(cell.optimal_score(), std::move(coordinate)));

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
    decltype(auto) track_last_row_cell(cell_t && cell, matrix_coordinate_type coordinate) noexcept
    {
        if (test_last_row_cell && !test_every_cell)
            binary_update_operation(std::forward_as_tuple(optimal_score, optimal_coordinate),
                                    std::forward_as_tuple(cell.optimal_score(), std::move(coordinate)));

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
    decltype(auto) track_last_column_cell(cell_t && cell, matrix_coordinate_type coordinate) noexcept
    {
        if (test_last_column_cell && !test_every_cell)
            binary_update_operation(std::forward_as_tuple(optimal_score, optimal_coordinate),
                                    std::forward_as_tuple(cell.optimal_score(), std::move(coordinate)));

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
    decltype(auto) track_final_cell(cell_t && cell, matrix_coordinate_type coordinate) noexcept
    {
        if (!(test_every_cell || test_last_row_cell || test_last_column_cell))
            binary_update_operation(std::forward_as_tuple(optimal_score, optimal_coordinate),
                                    std::forward_as_tuple(cell.optimal_score(), std::move(coordinate)));

        return std::forward<cell_t>(cell);
    }

    //!\brief Resets the optimum such that a new alignment can be computed.
    void reset_optimum() noexcept
    {
        optimal_score = std::numeric_limits<score_type>::lowest();
        optimal_coordinate = {};
    }
};
} // namespace seqan3::detail
