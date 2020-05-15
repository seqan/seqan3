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
#include <seqan3/std/ranges>

#include <seqan3/alignment/pairwise/detail/policy_optimum_tracker.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/type_traits/lazy.hpp>
#include <seqan3/range/views/zip.hpp>

namespace seqan3::detail
{

/*!\brief A binary operation to update the alignment optimum for the vectorised global alignment algorithm.
 * \ingroup pairwise_alignment
 *
 * \tparam simd_t The simd index type; must model seqan3::simd::simd_concept.
 *
 * \details
 *
 * Updates the current alignment optimum with the new score and the respective coordinate only if the coordinate
 * matches the stored column and row index for the respective alignment. The respective coordinates need to be
 * pre-computed.
 */
template <simd_concept simd_index_t>
struct alignment_optimum_updater_greater_equal_global_alignment_simd
{
    //!\brief The column indices that need to match.
    simd_index_t column_index{};
    //!\brief The row indices that need to match.
    simd_index_t row_index{};

    /*!\brief The binary compare and update operation.
     * \tparam lhs_t The type of the left hand side; must model seqan3::tuple_like with a tuple size of 2.
     * \tparam rhs_t The type of the right hand side; must model seqan3::tuple_like with a tuple size of 2.
     *
     * \param[in,out] optimal_score_coordinate_pair The current optimum.
     * \param[in] current_cell_score_coordinate_pair The new score and coordinate to compare with.
     *
     * \details
     *
     * Requires that first value of the tuple represents the simd score and the second type the matrix coordinate over
     * the simd column and row index.
     */
    template <tuple_like lhs_t, tuple_like rhs_t>
    //!\cond
        requires std::tuple_size_v<lhs_t> == 2 &&
                 std::tuple_size_v<rhs_t> == 2 &&
                 simd_concept<std::tuple_element_t<0, lhs_t>> &&
                 simd_concept<std::tuple_element_t<0, rhs_t>>
    //!\endcond
    void operator()(lhs_t && optimal_score_coordinate_pair,
                    rhs_t && current_cell_score_coordinate_pair) const
    {
        auto && [optimal_score, optimal_coordinate] = optimal_score_coordinate_pair;
        auto && [current_cell_score, current_cell_coordinate] = current_cell_score_coordinate_pair;

        auto mask = (column_index == current_cell_coordinate.col) && (row_index == current_cell_coordinate.row);
        optimal_score = (mask) ? current_cell_score : optimal_score;
        optimal_coordinate.col = (mask) ? column_index : optimal_coordinate.col;
        optimal_coordinate.row = (mask) ? row_index : optimal_coordinate.row;
    }
};

/*!\brief Implements the tracker to store the global optimum for a particular alignment computation.
 * \ingroup pairwise_alignment
 * \copydetails seqan3::detail::policy_optimum_tracker
 */
template <typename alignment_configuration_t, std::semiregular binary_update_operation_t>
//!\cond
    requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
//!\endcond
class policy_optimum_tracker_simd :
    protected policy_optimum_tracker<alignment_configuration_t, binary_update_operation_t>
{
protected:
    //!\brief The type of the base class.
    using base_policy_t = policy_optimum_tracker<alignment_configuration_t, binary_update_operation_t>;

    // Import the configured score type.
    using typename base_policy_t::traits_type;
    using typename base_policy_t::score_type;

    //!\brief The scalar type of the simd vector.
    using scalar_type = typename simd::simd_traits<score_type>::scalar_type;
    //!\brief The original non-simd score type.
    using original_score_type = typename traits_type::original_score_type;

    static_assert(simd_concept<score_type>, "Must be a simd type!");

    // Import base variables into class scope.
    using base_policy_t::binary_update_operation;
    using base_policy_t::optimal_score;
    using base_policy_t::optimal_coordinate;
    //!\brief The individual offsets used for padding the sequences.
    std::array<original_score_type, simd_traits<score_type>::length> padding_offsets{};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_optimum_tracker_simd() = default; //!< Defaulted.
    policy_optimum_tracker_simd(policy_optimum_tracker_simd const &) = default; //!< Defaulted.
    policy_optimum_tracker_simd(policy_optimum_tracker_simd &&) = default; //!< Defaulted.
    policy_optimum_tracker_simd & operator=(policy_optimum_tracker_simd const &) = default; //!< Defaulted.
    policy_optimum_tracker_simd & operator=(policy_optimum_tracker_simd &&) = default; //!< Defaulted.
    ~policy_optimum_tracker_simd() = default; //!< Defaulted.

    /*!\brief Construction and initialisation using the alignment configuration.
     * \param[in] config The alignment configuration (not used in this context).
     *
     * \details
     *
     * Resets the optimum on construction. Sets track_last_row_cell and track_last_column_cell to true in order
     * to check if this is a cell with a potential optimum.
     */
    policy_optimum_tracker_simd(alignment_configuration_t const & SEQAN3_DOXYGEN_ONLY(config))
    {
        base_policy_t::test_last_row_cell = true;
        base_policy_t::test_last_column_cell = true;

        reset_optimum();
    }
    //!\}

    //!\copydoc seqan3::detail::policy_optimum_tracker::reset_optimum
    void reset_optimum()
    {
        using index_t = typename traits_type::matrix_index_type;

        optimal_score = simd::fill<score_type>(std::numeric_limits<scalar_type>::lowest());
        optimal_coordinate.row = index_t{};
        optimal_coordinate.col = index_t{};
    }

    /*!\brief Initialises the tracker and possibly the binary update operation.
     * \tparam sequence1_collection_t The type of the sequence collection; must model std::ranges::input_range.
     * \tparam sequence2_collection_t The type of the sequence collection; must model std::ranges::input_range.
     *
     * \param[in] sequence1_collection The collection over sequences used for the initialisation of the tracker.
     * \param[in] sequence2_collection The collection over sequences used for the initialisation of the tracker.
     *
     * \details
     *
     * Initialises the binary update operation by pre-computing the coordinates for each individual matrix at which
     * the global alignment score can be found. Inside of the inter-sequence vectorisation layout the sequences might
     * have different sizes. Thus, the end coordinate for each individual alignment matrix can differ. The score is
     * populated to the end of the matrix and only the respective coordinates are tracked for the optimal score.
     * Finally, the added offset is removed from the score to obtain the true value.
     */
    template <std::ranges::input_range sequence1_collection_t, std::ranges::input_range sequence2_collection_t>
    void initialise_tracker(sequence1_collection_t & sequence1_collection,
                            sequence2_collection_t & sequence2_collection)
    {
        using index_t = typename traits_type::matrix_index_type;
        using scalar_index_t = typename simd_traits<index_t>::scalar_type;
        //
        scalar_index_t largest_sequence1_size{};
        scalar_index_t largest_sequence2_size{};
        alignas(alignof(index_t)) std::array<scalar_index_t, traits_type::alignments_per_vector> sequence1_sizes{};
        alignas(alignof(index_t)) std::array<scalar_index_t, traits_type::alignments_per_vector> sequence2_sizes{};

        // First, get all dimensions from the sequences and keep track of the maximal sizes.
        size_t index{};
        for (auto && [sequence1, sequence2] : views::zip(sequence1_collection, sequence2_collection))
        {
            sequence1_sizes[index] = static_cast<scalar_index_t>(std::ranges::distance(sequence1));
            sequence2_sizes[index] = static_cast<scalar_index_t>(std::ranges::distance(sequence2));
            largest_sequence1_size = std::max(largest_sequence1_size, sequence1_sizes[index]);
            largest_sequence2_size = std::max(largest_sequence2_size, sequence2_sizes[index]);
            ++index;
        }

        assert(index > 0); // We do not expect empty sequence collections here.

        // Second, get the offset for each individual end coordinate to project the cell to the last row or
        // column of the global alignment matrix. Choose the smallest distance, since this gives the correct offset
        // to the projected end cell.
        do
        {
            --index;

            assert(sequence1_sizes[index] <= largest_sequence1_size);
            assert(sequence2_sizes[index] <= largest_sequence2_size);

            padding_offsets[index] = std::min((largest_sequence1_size - sequence1_sizes[index]),
                                                (largest_sequence2_size - sequence2_sizes[index]));
            sequence1_sizes[index] += padding_offsets[index];
            sequence2_sizes[index] += padding_offsets[index];
        } while (index > 0);

        // Load the target coordinate indices from the respective arrays.
        binary_update_operation.column_index = simd::load<index_t>(sequence1_sizes.data());
        binary_update_operation.row_index = simd::load<index_t>(sequence2_sizes.data());
    }
};
} // namespace seqan3::detail
