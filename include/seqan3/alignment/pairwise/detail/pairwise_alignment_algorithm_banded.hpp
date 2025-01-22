// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::pairwise_alignment_algorithm.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <ranges>

#include <seqan3/alignment/pairwise/detail/pairwise_alignment_algorithm.hpp>
#include <seqan3/utility/views/slice.hpp>

namespace seqan3::detail
{

/*!\brief The alignment algorithm type to compute the banded standard pairwise alignment using dynamic programming.
 * \implements std::invocable
 * \ingroup alignment_pairwise
 * \copydetails seqan3::detail::pairwise_alignment_algorithm
 */
template <typename alignment_configuration_t, typename... policies_t>
    requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
class pairwise_alignment_algorithm_banded :
    protected pairwise_alignment_algorithm<alignment_configuration_t, policies_t...>
{
protected:
    //!\brief The alignment configuration traits type with auxiliary information extracted from the configuration type.
    using base_algorithm_t = pairwise_alignment_algorithm<alignment_configuration_t, policies_t...>;

    // Import types from base class.
    using typename base_algorithm_t::alignment_result_type;
    using typename base_algorithm_t::score_type;
    using typename base_algorithm_t::traits_type;

    static_assert(!std::same_as<alignment_result_type, empty_type>, "Alignment result type was not configured.");
    static_assert(traits_type::is_banded, "Alignment configuration must have band configured.");

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    pairwise_alignment_algorithm_banded() = default;                                            //!< Defaulted.
    pairwise_alignment_algorithm_banded(pairwise_alignment_algorithm_banded const &) = default; //!< Defaulted.
    pairwise_alignment_algorithm_banded(pairwise_alignment_algorithm_banded &&) = default;      //!< Defaulted.
    pairwise_alignment_algorithm_banded &
    operator=(pairwise_alignment_algorithm_banded const &) = default;                                  //!< Defaulted.
    pairwise_alignment_algorithm_banded & operator=(pairwise_alignment_algorithm_banded &&) = default; //!< Defaulted.
    ~pairwise_alignment_algorithm_banded() = default;                                                  //!< Defaulted.

    /*!\brief Constructs and initialises the algorithm using the alignment configuration.
     * \param config The configuration passed into the algorithm.
     *
     * \details
     *
     * Initialises the algorithm given the user settings from the alignment configuration object.
     *
     * \throws seqan3::invalid_alignment_configuration.
     */
    pairwise_alignment_algorithm_banded(alignment_configuration_t const & config) : base_algorithm_t(config)
    {}
    //!\}

    /*!\name Invocation
     * \{
     */
    //!\copydoc seqan3::detail::pairwise_alignment_algorithm::operator()(indexed_sequence_pairs_t && indexed_sequence_pairs, callback_t && callback)
    template <indexed_sequence_pair_range indexed_sequence_pairs_t, typename callback_t>
        requires std::invocable<callback_t, alignment_result_type>
    void operator()(indexed_sequence_pairs_t && indexed_sequence_pairs, callback_t && callback)
    {
        using std::get;

        for (auto && [sequence_pair, idx] : indexed_sequence_pairs)
        {
            size_t sequence1_size = std::ranges::distance(get<0>(sequence_pair));
            size_t const sequence2_size = std::ranges::distance(get<1>(sequence_pair));

            auto && [alignment_matrix, index_matrix] =
                this->acquire_matrices(sequence1_size, sequence2_size, this->lowest_viable_score());

            // Initialise the cell updater with the dimensions of the regular matrix.
            this->compare_and_set_optimum.set_target_indices(row_index_type{sequence2_size},
                                                             column_index_type{sequence1_size});

            // Shrink the first sequence if the band ends before its actual end.
            sequence1_size = std::min(sequence1_size, this->upper_diagonal + sequence2_size);

            using sequence1_difference_t = std::ranges::range_difference_t<decltype(get<0>(sequence_pair))>;

            compute_matrix(std::views::take(get<0>(sequence_pair), static_cast<sequence1_difference_t>(sequence1_size)),
                           get<1>(sequence_pair),
                           alignment_matrix,
                           index_matrix);
            this->make_result_and_invoke(std::forward<decltype(sequence_pair)>(sequence_pair),
                                         std::move(idx),
                                         this->optimal_score,
                                         this->optimal_coordinate,
                                         alignment_matrix,
                                         callback);
        }
    }

    //!\overload
    template <indexed_sequence_pair_range indexed_sequence_pairs_t, typename callback_t>
        requires traits_type::is_vectorised && std::invocable<callback_t, alignment_result_type>
    auto operator()(indexed_sequence_pairs_t && indexed_sequence_pairs, callback_t && callback)
    {
        using simd_collection_t = std::vector<score_type, aligned_allocator<score_type, alignof(score_type)>>;
        using original_score_t = typename traits_type::original_score_type;

        // Extract the batch of sequences for the first and the second sequence.
        auto seq1_collection = indexed_sequence_pairs | views::elements<0> | views::elements<0>;
        auto seq2_collection = indexed_sequence_pairs | views::elements<0> | views::elements<1>;

        this->initialise_tracker(seq1_collection, seq2_collection);

        // Convert batch of sequences to sequence of simd vectors.
        thread_local simd_collection_t simd_seq1_collection{};
        thread_local simd_collection_t simd_seq2_collection{};

        this->convert_batch_of_sequences_to_simd_vector(simd_seq1_collection,
                                                        seq1_collection,
                                                        this->scoring_scheme.padding_symbol);
        this->convert_batch_of_sequences_to_simd_vector(simd_seq2_collection,
                                                        seq2_collection,
                                                        this->scoring_scheme.padding_symbol);

        size_t const sequence1_size = std::ranges::distance(simd_seq1_collection);
        size_t const sequence2_size = std::ranges::distance(simd_seq2_collection);

        auto && [alignment_matrix, index_matrix] =
            this->acquire_matrices(sequence1_size, sequence2_size, this->lowest_viable_score());

        compute_matrix(simd_seq1_collection, simd_seq2_collection, alignment_matrix, index_matrix);

        size_t index = 0;
        for (auto && [sequence_pair, idx] : indexed_sequence_pairs)
        {
            original_score_t score = this->optimal_score[index]
                                   - (this->padding_offsets[index] * this->scoring_scheme.padding_match_score());
            matrix_coordinate coordinate{row_index_type{size_t{this->optimal_coordinate.row[index]}},
                                         column_index_type{size_t{this->optimal_coordinate.col[index]}}};
            this->make_result_and_invoke(std::forward<decltype(sequence_pair)>(sequence_pair),
                                         std::move(idx),
                                         std::move(score),
                                         std::move(coordinate),
                                         alignment_matrix,
                                         callback);
            ++index;
        }
    }
    //!\}

protected:
    /*!\brief Compute the actual banded alignment.
     * \tparam sequence1_t The type of the first sequence; must model std::ranges::forward_range.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::forward_range.
     * \tparam alignment_matrix_t The type of the alignment matrix; must model std::ranges::input_range and its
     *                            std::ranges::range_reference_t type must model std::ranges::forward_range.
     * \tparam index_matrix_t The type of the index matrix; must model std::ranges::input_range and its
     *                            std::ranges::range_reference_t type must model std::ranges::forward_range.
     *
     * \param[in] sequence1 The first sequence to compute the alignment for.
     * \param[in] sequence2 The second sequence to compute the alignment for.
     * \param[in] alignment_matrix The alignment matrix to compute.
     * \param[in] index_matrix The index matrix corresponding to the alignment matrix.
     *
     * \details
     *
     * In the banded alignment the iteration of the inner columns is split into two phases. The first phase
     * reuses the unbanded column computation and assumes that the banded score matrix always starts at the beginning
     * of the matrix. In the second phase the special interface
     * seqan3::detail::pairwise_alignment_algorithm_banded::compute_band_column is used to compute the banded column.
     * The current implementation (this might change when we need to work with full-matrices like Waterman-Eggert does)
     * assumes that the first cell of the current score matrix column is always the first
     * cell to compute in every column. The respective seqan3::detail::score_matrix_single_column is resized to the band
     * size plus one additional field to cover the end of the band where the end of the column is not reached yet.
     * This cell will never be written to but only read from, i.e. it represents minus infinity. This allows the
     * algorithm to reuse the standard unbanded cell computation. The following figure depicts the referenced
     * cell of the underlying score matrix (assuming seqan3::detail::score_matrix_single_column):
     *
     * ```
     *       A G G T C A
     *     0 1 2 3 4 5 6
     *    |–|—|—|—|—|—|—|
     *   0|0|0|0|0|0| | |
     * A 1|1|1|1|1|1|0| |
     * C 2|x|2|2|2|2|1|0|
     * G 3| |x|3|3|3|2|1|
     * T 4| | |x|4|4|3|2|
     *```

     * The coordinate matrix represents the global matrix index and not the local band coordinate. Data structures that
     * require the coordinate might need to map the global matrix coordinate to their local coordinate:
     *
     * ```
     *             A     G     G     T     C     A
     *       0     1     2     3     4     5     6
     *    |–––––|–––––|–––––|–––––|–––––|–––––|–––––|
     *   0|(0,0)|(0,1)|(0,2)|(0,3)|(0,4)|     |     |
     * A 1|(1,0)|(1,1)|(1,2)|(1,3)|(1,4)|(1,5)|     |
     * C 2|     |(2,1)|(2,2)|(2,3)|(2,4)|(2,5)|(2,6)|
     * G 3|     |     |(3,2)|(3,3)|(3,4)|(3,5)|(3,6)|
     * T 4|     |     |     |(4,3)|(4,4)|(4,5)|(4,6)|
     *```
     */
    template <std::ranges::forward_range sequence1_t,
              std::ranges::forward_range sequence2_t,
              std::ranges::input_range alignment_matrix_t,
              std::ranges::input_range index_matrix_t>
        requires std::ranges::forward_range<std::ranges::range_reference_t<alignment_matrix_t>>
              && std::ranges::forward_range<std::ranges::range_reference_t<index_matrix_t>>
    void compute_matrix(sequence1_t && sequence1,
                        sequence2_t && sequence2,
                        alignment_matrix_t && alignment_matrix,
                        index_matrix_t && index_matrix)
    {
        // ---------------------------------------------------------------------
        // Initialisation phase: allocate memory and initialise first column.
        // ---------------------------------------------------------------------

        this->reset_optimum(); // Reset the tracker for the new alignment computation.

        auto alignment_matrix_it = alignment_matrix.begin();
        auto indexed_matrix_it = index_matrix.begin();

        using row_index_t = std::ranges::range_difference_t<sequence2_t>;    // row_size = |sequence2| + 1
        using column_index_t = std::ranges::range_difference_t<sequence1_t>; // column_size = |sequence1| + 1

        row_index_t row_size = std::max<int32_t>(0, -this->lower_diagonal);
        column_index_t const column_size = std::max<int32_t>(0, this->upper_diagonal);
        this->initialise_column(*alignment_matrix_it, *indexed_matrix_it, std::views::take(sequence2, row_size));

        // ---------------------------------------------------------------------
        // 1st recursion phase: band intersects with the first row.
        // ---------------------------------------------------------------------

        for (auto alphabet1 : std::views::take(sequence1, column_size))
        {
            this->compute_column(*++alignment_matrix_it,
                                 *++indexed_matrix_it,
                                 alphabet1,
                                 std::views::take(sequence2, ++row_size));
        }

        // ---------------------------------------------------------------------
        // 2nd recursion phase: iterate until the end of the matrix.
        // ---------------------------------------------------------------------

        row_index_t first_row_index = 0u;
        for (auto alphabet1 : std::views::drop(sequence1, column_size))
        {
            compute_band_column(*++alignment_matrix_it,
                                std::views::drop(*++indexed_matrix_it, first_row_index + 1),
                                alphabet1,
                                views::slice(sequence2, first_row_index, ++row_size));
            ++first_row_index;
        }

        // ---------------------------------------------------------------------
        // Final phase: track score of last column
        // ---------------------------------------------------------------------

        auto alignment_column = *alignment_matrix_it;
        auto cell_index_column = std::views::drop(*indexed_matrix_it, first_row_index);

        auto alignment_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();

        this->track_last_column_cell(*alignment_column_it, *cell_index_column_it);

        for (row_index_t last_row = std::min<row_index_t>(std::ranges::distance(sequence2), row_size);
             first_row_index < last_row;
             ++first_row_index)
            this->track_last_column_cell(*++alignment_column_it, *++cell_index_column_it);

        this->track_final_cell(*alignment_column_it, *cell_index_column_it);
    }

    /*!\brief Computes a column of the band that does not start in the first row of the alignment matrix.
     * \tparam alignment_column_t The type of the alignment column; must model std::ranges::forward_range.
     * \tparam cell_index_column_t The type of the indexed column; must model std::ranges::input_range.
     * \tparam alphabet1_t The type of the current symbol of sequence1.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::input_range.
     *
     * \param[in] alignment_column The current alignment matrix column to compute.
     * \param[in] cell_index_column The current index matrix column to get the respective cell indices.
     * \param[in] alphabet1 The current symbol of sequence1.
     * \param[in] sequence2 The second sequence to align against `alphabet1`.
     *
     * \details
     *
     * Computes the alignment for the given alignment matrix column. The function splits the computation of the column
     * into three phases: the initialisation phase, the iteration phase, and the final phase. In the initialisation
     * phase the first cell of the column is computed and in the iteration phase all remaining cells are computed.
     * In the final phase the last cell is possibly evaluated for a new alignment optimum.
     * Note that the length of `sequence2` determines the size of the column.
     *
     * ### Implementation with a score matrix using linear memory
     *
     * When the score matrix uses only linear memory, i.e. only one column is stored, the algorithm reuses the cells of
     * the same column to compute the current one. This means, before the current cell is updated its value is
     * cached and used as the previous diagonal value when computing the next cell below.
     * In the banded case, however, the position of the referenced cell is shifted by one and needs a different
     * handling.
     *
     * ```
     *    0 1 2 3 4 5 6
     *    _____________
     * 0 |0|0|0|\      |
     * 1 |1|1|1|0|\    |
     * 2 |\|2|2|1|0|\  |
     * 3 |  \|3|2|1|0|\|
     * 4 |    \|3|2|1|0|
     * 5 |      \|3|2|1|
     *    –––––––––––––
     *
     * ```
     * The picture above depicts the banded matrix with the indices of the column. As long as the band touches
     * the first row (column 0 - 2), the indices of the actual stored column refer to the same position as the
     * previous column. Hence, to compute these columns the regular algorithm
     * seqan3::detail::pairwise_alignment_algorithm::compute_column can be used.
     * Starting at column 3 the first cell of the band moves downwards, causing a shift in the position of the previous
     * cell. In this state, the value of the current cell represents the previous diagonal value before it is updated.
     * To read the previous horizontal value the next cell below has to be dereferenced.
     * Accordingly, two iterators are used to point to the respective cells in the matrix. The first one points to the
     * current cell (the one that is written to) and the second points to the next cell (the one where the
     * horizontal and vertical scores are read from). After computing the last cell of the column the value of the
     * current iterator can be used to track the score of the cell.
     */
    template <std::ranges::forward_range alignment_column_t,
              std::ranges::input_range cell_index_column_t,
              typename alphabet1_t,
              std::ranges::input_range sequence2_t>
    void compute_band_column(alignment_column_t && alignment_column,
                             cell_index_column_t && cell_index_column,
                             alphabet1_t const & alphabet1,
                             sequence2_t && sequence2)
    {
        // ---------------------------------------------------------------------
        // Initial phase: prepare column and initialise first cell
        // ---------------------------------------------------------------------

        auto current_alignment_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();

        // Points to the last valid cell in the column.
        decltype(current_alignment_column_it) next_alignment_column_it{current_alignment_column_it};
        auto cell = *current_alignment_column_it;
        cell = this->track_cell(
            this->initialise_band_first_cell(cell.best_score(),
                                             *++next_alignment_column_it,
                                             this->scoring_scheme.score(alphabet1, *std::ranges::begin(sequence2))),
            *cell_index_column_it);

        // ---------------------------------------------------------------------
        // Iteration phase: iterate over column and compute each cell
        // ---------------------------------------------------------------------

        for (auto && alphabet2 : std::views::drop(sequence2, 1))
        {
            current_alignment_column_it = next_alignment_column_it;
            auto cell = *current_alignment_column_it;
            cell = this->track_cell(this->compute_inner_cell(cell.best_score(),
                                                             *++next_alignment_column_it,
                                                             this->scoring_scheme.score(alphabet1, alphabet2)),
                                    *++cell_index_column_it);
        }

        // ---------------------------------------------------------------------
        // Final phase: track last cell
        // ---------------------------------------------------------------------

        this->track_last_row_cell(*current_alignment_column_it, *cell_index_column_it);
    }
};

} // namespace seqan3::detail
