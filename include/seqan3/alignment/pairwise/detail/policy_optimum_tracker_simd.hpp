// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::policy_optimum_tracker_simd.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <limits>
#include <ranges>

#include <seqan3/alignment/pairwise/detail/policy_optimum_tracker.hpp>
#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/simd.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>
#include <seqan3/utility/type_traits/lazy_conditional.hpp>
#include <seqan3/utility/views/zip.hpp>

namespace seqan3::detail
{

/*!\brief Function object that compares and updates the alignment optimum for the vectorised global alignment algorithm.
 * \ingroup alignment_pairwise
 *
 * \details
 *
 * This operation is specific to the global alignment in vectorised mode. The tracking of the last cells of the
 * different alignment matrices that are computed simultaneously in one vector unit depends on how the scoring of the
 * global alignment works. Any alignment matrix that is smaller than the largest matrix defined by the longest sequence
 * in the collection of the first sequences and in the collection of the second sequences, will have its last cell
 * in the middle of the encompassing matrix. In order to track this cell without checking every cell of the alignment
 * matrix for the correct coordinate, the last cell of every contained matrix is projected along its diagonal to either
 * the last row or the last column of the encompassing matrix. Within the algorithm only the cells of the last row,
 * respectively column are tracked. Since the cells to track are fixed, the respective coordinates for every
 * contained matrix can be precomputed. Subsequently, within this update operation the score is only updated if the
 * coordinate of the current cell compares equal to the precomputed coordinate of any of the contained matrices.
 */
struct max_score_updater_simd_global
{

    /*!\brief Compares and updates the optimal score-coordinate pair.
     * \tparam score_t The type of the score to track; must model std::assignable_from `const & score_t`.
     * \tparam coordinate_t The type of the coordinate to track; must be a seqan3::matrix_index type with members that
     *                      model seqan3::simd::simd_concept.
     *
     * \param[in,out] optimal_score The optimal score to update.
     * \param[in,out] optimal_coordinate The optimal coordinate to update.
     * \param[in] current_score The score of the current cell.
     * \param[in] current_coordinate The coordinate of the current cell.
     *
     * \details
     *
     * Compares the coordinate of the current cell with the precomputed coordinate which represent the projected
     * last cells of the contained matrices. If a coordinate matches the precomputed one, then the respective score
     * will be set for the optimal score.
     */
    template <typename score_t, typename coordinate_t>
        requires (std::assignable_from<score_t &, score_t const &> &&
                  requires (coordinate_t coordinate)
                  {
                      requires simd_concept<decltype(coordinate.col)>;
                      requires simd_concept<decltype(coordinate.row)>;
                  })
    void operator()(score_t & optimal_score,
                    coordinate_t const & optimal_coordinate,
                    score_t current_score,
                    coordinate_t const & current_coordinate) const noexcept
    {
        auto mask =
            (optimal_coordinate.col == current_coordinate.col) && (optimal_coordinate.row == current_coordinate.row);
        optimal_score = (mask) ? std::move(current_score) : optimal_score;
    }
};

/*!\brief Implements the tracker to store the global optimum for a particular alignment computation.
 * \ingroup alignment_pairwise
 * \copydetails seqan3::detail::policy_optimum_tracker
 */
template <typename alignment_configuration_t, std::semiregular optimum_updater_t>
    requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
          && std::invocable<
                 optimum_updater_t,
                 typename alignment_configuration_traits<alignment_configuration_t>::score_type &,
                 typename alignment_configuration_traits<alignment_configuration_t>::matrix_coordinate_type &,
                 typename alignment_configuration_traits<alignment_configuration_t>::score_type,
                 typename alignment_configuration_traits<alignment_configuration_t>::matrix_coordinate_type>
class policy_optimum_tracker_simd : protected policy_optimum_tracker<alignment_configuration_t, optimum_updater_t>
{
protected:
    //!\brief The type of the base class.
    using base_policy_t = policy_optimum_tracker<alignment_configuration_t, optimum_updater_t>;

    // Import the configured score type.
    using typename base_policy_t::score_type;
    using typename base_policy_t::traits_type;

    //!\brief The scalar type of the simd vector.
    using scalar_type = typename simd::simd_traits<score_type>::scalar_type;
    //!\brief The original non-simd score type.
    using original_score_type = typename traits_type::original_score_type;

    static_assert(simd_concept<score_type>, "Must be a simd type!");

    // Import base variables into class scope.
    using base_policy_t::compare_and_set_optimum;
    using base_policy_t::optimal_coordinate;
    using base_policy_t::optimal_score;
    //!\brief The individual offsets used for padding the sequences.
    std::array<original_score_type, simd_traits<score_type>::length> padding_offsets{};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_optimum_tracker_simd() = default;                                                //!< Defaulted.
    policy_optimum_tracker_simd(policy_optimum_tracker_simd const &) = default;             //!< Defaulted.
    policy_optimum_tracker_simd(policy_optimum_tracker_simd &&) = default;                  //!< Defaulted.
    policy_optimum_tracker_simd & operator=(policy_optimum_tracker_simd const &) = default; //!< Defaulted.
    policy_optimum_tracker_simd & operator=(policy_optimum_tracker_simd &&) = default;      //!< Defaulted.
    ~policy_optimum_tracker_simd() = default;                                               //!< Defaulted.

    /*!\brief Construction and initialisation using the alignment configuration.
     * \param[in] config The alignment configuration (not used in this context).
     *
     * \details
     *
     * Initialises the object to always track the last row and column, since this is needed for the vectorised global
     * alignment.
     */
    policy_optimum_tracker_simd(alignment_configuration_t const & config) : base_policy_t{config}
    {
        base_policy_t::test_last_row_cell = true;
        base_policy_t::test_last_column_cell = true;
    }
    //!\}

    //!\copydoc seqan3::detail::policy_optimum_tracker::reset_optimum
    void reset_optimum()
    {
        optimal_score = simd::fill<score_type>(std::numeric_limits<scalar_type>::lowest());
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
     * Initialises the binary max score operation by pre-computing the coordinates for each individual matrix at which
     * the global alignment score can be found. Inside of the inter-sequence vectorisation layout the sequences might
     * have different sizes. Thus, the end coordinate for each individual alignment matrix can differ. The score is
     * populated to the end of the matrix and only the respective coordinates are tracked for the optimal score.
     * Finally, the added offset is removed from the score to obtain the true value.
     *
     * ### Example
     *
     * Consider the following collections:
     * collection 1: ["aaa", "aa",   "a"]
     * collection 2: [  "a", "aa", "aaa"]
     *
     * Based on the length of the sequences the encompassing alignment matrix has the dimensions 4x4.
     * The following graphic depicts this matrix. The number `i` marks the last cell for the contained matrix given the
     * sequence pair from above: 1 -> (3,1); 2 -> (2,2); 3 -> (1,3)
     * ```
     *   |0|1|2|3|
     *  -|-|-|-|-|
     *  0| | | | |
     *  -|-|-|-|-|
     *  1| | | |1|
     *  -|-|-|-|-|
     *  2| | |2| |
     *  -|-|-|-|-|
     *  3| |3| | |
     * ```
     *
     * As can be seen, the end of each matrix must not necessarily be the end of the encompassing matrix. To avoid
     * tracking every cell in the matrix the end point candidates will be projected along the diagonal to the
     * last row or column of the encompassing matrix. In the example above, 1 and 3 are projected with an offset of `0`,
     * while 2 is projected to the coordinate (3,3) with an projection offset of `1`. During the computation of the
     * alignment the simd scoring scheme ensures that outside of the original matrix a fixed cost is added to every
     * cell. The original score is shifted by this cost multiplied with the computed projection offset.
     * In the global alignment it is suffcient to only track the optimal score in the last row and column of the
     * encompassing matrix and only at the precomputed coordinate projections. Eventually, the score offset is
     * subtracted to obtain the original score.
     */
    template <std::ranges::input_range sequence1_collection_t, std::ranges::input_range sequence2_collection_t>
    void initialise_tracker(sequence1_collection_t & sequence1_collection,
                            sequence2_collection_t & sequence2_collection)
    {
        using index_t = typename traits_type::matrix_index_type;
        using scalar_index_t = typename simd_traits<index_t>::scalar_type;

        scalar_index_t largest_sequence1_size{};
        scalar_index_t largest_sequence2_size{};
        alignas(alignof(index_t)) std::array<scalar_index_t, traits_type::alignments_per_vector> sequence1_sizes{};
        alignas(alignof(index_t)) std::array<scalar_index_t, traits_type::alignments_per_vector> sequence2_sizes{};

        // First, get all dimensions from the sequences and keep track of the maximal size in either dimension.
        size_t sequence_count{};
        for (auto && [sequence1, sequence2] : views::zip(sequence1_collection, sequence2_collection))
        {
            sequence1_sizes[sequence_count] = std::ranges::distance(sequence1);
            sequence2_sizes[sequence_count] = std::ranges::distance(sequence2);
            largest_sequence1_size = std::max(largest_sequence1_size, sequence1_sizes[sequence_count]);
            largest_sequence2_size = std::max(largest_sequence2_size, sequence2_sizes[sequence_count]);
            ++sequence_count;
        }

        // Second, determine the offset for each individual end-coordinate which is used to project the cell to the
        // last row or column of the global alignment matrix. Choose the smallest distance as the correct offset
        // to the projected cell.
        for (size_t index = 0; index != sequence_count; ++index)
        {
            assert(sequence1_sizes[index] <= largest_sequence1_size);
            assert(sequence2_sizes[index] <= largest_sequence2_size);

            padding_offsets[index] = std::min(largest_sequence1_size - sequence1_sizes[index],
                                              largest_sequence2_size - sequence2_sizes[index]);
            sequence1_sizes[index] += padding_offsets[index];
            sequence2_sizes[index] += padding_offsets[index];
        }

        // Load the target coordinate indices from the respective arrays.
        optimal_coordinate.col = simd::load<index_t>(sequence1_sizes.data());
        optimal_coordinate.row = simd::load<index_t>(sequence2_sizes.data());
    }
};
} // namespace seqan3::detail
