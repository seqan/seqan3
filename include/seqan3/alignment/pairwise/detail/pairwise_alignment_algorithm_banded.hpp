// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::pairwise_alignment_algorithm.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

#include <seqan3/alignment/exception.hpp>
#include <seqan3/alignment/pairwise/detail/pairwise_alignment_algorithm.hpp>
#include <seqan3/range/views/drop.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/range/views/slice.hpp>

namespace seqan3::detail
{

/*!\brief The alignment algorithm type to compute the banded standard pairwise alignment using dynamic programming.
 * \implements std::invocable
 * \ingroup pairwise_alignment
 * \copydetails seqan3::detail::pairwise_alignment_algorithm
 */
template <typename alignment_configuration_t, typename ...policies_t>
//!\cond
    requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
//!\endcond
class pairwise_alignment_algorithm_banded :
    protected pairwise_alignment_algorithm<alignment_configuration_t, policies_t...>
{
protected:
    //!\brief The alignment configuration traits type with auxiliary information extracted from the configuration type.
    using base_algorithm_t = pairwise_alignment_algorithm<alignment_configuration_t, policies_t...>;

    // import types from base class.
    using typename base_algorithm_t::traits_type;
    using typename base_algorithm_t::scoring_scheme_type;
    using typename base_algorithm_t::alignment_result_type;

    static_assert(!std::same_as<alignment_result_type, empty_type>, "Alignment result type was not configured.");
    static_assert(traits_type::is_banded, "Alignment configuration must have band configured.");

    //!\brief The selected lower diagonal.
    int32_t lower_diagonal{};
    //!\brief The selected upper diagonal.
    int32_t upper_diagonal{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    pairwise_alignment_algorithm_banded() = default; //!< Defaulted.
    pairwise_alignment_algorithm_banded(pairwise_alignment_algorithm_banded const &) = default; //!< Defaulted.
    pairwise_alignment_algorithm_banded(pairwise_alignment_algorithm_banded &&) = default; //!< Defaulted.
    pairwise_alignment_algorithm_banded & operator=(pairwise_alignment_algorithm_banded const &) = default; //!< Defaulted.
    pairwise_alignment_algorithm_banded & operator=(pairwise_alignment_algorithm_banded &&) = default; //!< Defaulted.
    ~pairwise_alignment_algorithm_banded() = default; //!< Defaulted.

    /*!\brief Constructs and initialises the algorithm using the alignment configuration.
     * \param config The configuration passed into the algorithm.
     *
     * \details
     *
     * Initialises the algorithm given the user settings from the alignment configuration object.
     *
     * \throws seqan3::invalid_alignment_configuration if the given band settings are invalid.
     */
    pairwise_alignment_algorithm_banded(alignment_configuration_t const & config) : base_algorithm_t(config)
    {
        using seqan3::get;

        auto band = get<seqan3::align_cfg::band_fixed_size>(config);

        lower_diagonal = band.lower_diagonal.get();
        upper_diagonal = band.upper_diagonal.get();

        if (upper_diagonal < lower_diagonal || upper_diagonal < 0 || lower_diagonal > 0)
            throw invalid_alignment_configuration{"The selected band [" + std::to_string(lower_diagonal) + ":" +
                                                  std::to_string(upper_diagonal) + "] is not supported because it "
                                                  "does not completely encloses both sequences."};
    }
    //!\}

    /*!\name Invocation
     * \{
     */
    //!\copydoc seqan3::detail::pairwise_alignment_algorithm::operator()(indexed_sequence_pairs_t && indexed_sequence_pairs, callback_t && callback)
    template <indexed_sequence_pair_range indexed_sequence_pairs_t, typename callback_t>
    //!\cond
        requires std::invocable<callback_t, alignment_result_type>
    //!\endcond
    void operator()(indexed_sequence_pairs_t && indexed_sequence_pairs, callback_t && callback)
    {
        using result_value_t = typename alignment_result_value_type_accessor<alignment_result_type>::type;
        using std::get;

        for (auto && [sequence_pair, idx] : indexed_sequence_pairs)
        {
            result_value_t res{};
            res.id = idx;

            compute_matrix(get<0>(sequence_pair), get<1>(sequence_pair));
            res.score = this->tracked_optimum();
            callback(alignment_result_type{res});
        }
    }
    //!\}

protected:
    /*!\brief Checks whether the band is valid for the given sequence sizes.
     *
     * \param[in] sequence1_size The size of the first sequence.
     * \param[in] sequence2_size The size of the second sequence.
     *
     * \throws seqan3::invalid_alignment_configuration if the band is invalid for the given sequence sizes and the
     *         alignment configuration.
     */
    void check_valid_band_configuration(size_t const sequence1_size, size_t const sequence2_size) const
    {
        bool const upper_diagonal_ends_before_last_cell = (upper_diagonal + sequence2_size) < sequence1_size;
        bool const lower_diagonal_ends_behind_last_cell = (-lower_diagonal + sequence1_size) < sequence2_size;

        if (upper_diagonal_ends_before_last_cell || lower_diagonal_ends_behind_last_cell)
            throw invalid_alignment_configuration{"The selected band [" + std::to_string(lower_diagonal) + ":" +
                                                  std::to_string(upper_diagonal) + "] does not cover the last cell."};
    }

    /*!\brief Compute the actual banded alignment.
     * \tparam sequence1_t The type of the first sequence; must model std::ranges::forward_range.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::forward_range.
     *
     * \param[in] sequence1 The first sequence to compute the alignment for.
     * \param[in] sequence2 The second sequence to compute the alignment for.
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
     *
     * \throws seqan3::invalid_alignment_configuration if the band is not valid for the given sequences or
     *         std::bad_alloc if too much memory is allocated.
     */
    template <std::ranges::forward_range sequence1_t, std::ranges::forward_range sequence2_t>
    void compute_matrix(sequence1_t && sequence1, sequence2_t && sequence2)
    {
        size_t sequence1_size = std::ranges::distance(sequence1);
        size_t sequence2_size = std::ranges::distance(sequence2);

        check_valid_band_configuration(sequence1_size, sequence2_size);

        // ---------------------------------------------------------------------
        // Initialisation phase: allocate memory and initialise first column.
        // ---------------------------------------------------------------------

        this->reset_optimum(); // Reset the tracker for the new alignment computation.

        thread_local score_matrix_single_column<int32_t> local_score_matrix{};
        coordinate_matrix<uint32_t> local_index_matrix{};

        size_t const number_of_columns = sequence1_size + 1;
        size_t const number_of_rows = sequence2_size + 1;
        size_t const band_size = upper_diagonal - lower_diagonal + 1;

        local_score_matrix.resize(column_index_type{number_of_columns},
                                  row_index_type{band_size + 1},
                                  this->lowest_viable_score());
        local_index_matrix.resize(column_index_type{number_of_columns}, row_index_type{number_of_rows});

        auto alignment_matrix_it = local_score_matrix.begin();
        auto indexed_matrix_it = local_index_matrix.begin();

        size_t row_size = std::max<int32_t>(0, -lower_diagonal);
        size_t const column_size = std::max<int32_t>(0, upper_diagonal);
        this->initialise_column(*alignment_matrix_it, *indexed_matrix_it, sequence2 | views::take(row_size));

        // ---------------------------------------------------------------------
        // 1st recursion phase: band intersects with the first row.
        // ---------------------------------------------------------------------

        for (auto sequence1_value : sequence1 | views::take(column_size))
        {
            this->compute_column(*++alignment_matrix_it,
                                 *++indexed_matrix_it,
                                 sequence1_value,
                                 sequence2 | views::take(++row_size));
        }

        // ---------------------------------------------------------------------
        // 2nd recursion phase: iterate until the end of the matrix.
        // ---------------------------------------------------------------------

        size_t first_row_index = 0u;
        for (auto sequence1_value : sequence1 | views::drop(column_size))
        {
            compute_band_column(*++alignment_matrix_it,
                                *++indexed_matrix_it | views::drop(first_row_index + 1),
                                sequence1_value,
                                sequence2 | views::slice(first_row_index, ++row_size));
            ++first_row_index;
        }

        // ---------------------------------------------------------------------
        // Final phase: track score of last column
        // ---------------------------------------------------------------------

        auto alignment_column = *alignment_matrix_it;
        auto cell_index_column = *indexed_matrix_it | views::drop(first_row_index);

        auto alignment_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();

        this->track_last_column_cell(*alignment_column_it, *cell_index_column_it);

        for (; first_row_index < sequence2_size; ++first_row_index)
            this->track_last_column_cell(*++alignment_column_it, *++cell_index_column_it);

        this->track_final_cell(*alignment_column_it, *cell_index_column_it);
    }

    /*!\brief Computes a column of the band that does not start in the first row of the alignment matrix.
     * \tparam alignment_column_t The type of the alignment column; must model std::ranges::input_range.
     * \tparam cell_index_column_t The type of the indexed column; must model std::ranges::input_range.
     * \tparam sequence1_value_t The value type of sequence1; must model seqan3::semialphabet.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::input_range.
     *
     * \param[in] alignment_column The current alignment matrix column to compute.
     * \param[in] cell_index_column The current index matrix column to get the respective cell indices.
     * \param[in] sequence1_value The current symbol of sequence1.
     * \param[in] sequence2 The second sequence to align against `sequence1_value`.
     *
     * \details
     *
     * Computes the alignment for the given alignment matrix column. The function splits the computation of the column
     * into three phases: the initialisation phase, the iteration phase, and the final phase. In the initialisation
     * phase the first cell of the column is computed and in the iteration phase all remaining cells are computed.
     * In the final phase the last cell is possibly evaluated for a new alignment optimum.
     * Note that the length of `sequence2` determines the size of the column.
     */
    template <std::ranges::input_range alignment_column_t,
              std::ranges::input_range cell_index_column_t,
              semialphabet sequence1_value_t,
              std::ranges::input_range sequence2_t>
    void compute_band_column(alignment_column_t && alignment_column,
                             cell_index_column_t && cell_index_column,
                             sequence1_value_t const & sequence1_value,
                             sequence2_t && sequence2)
    {
        // ---------------------------------------------------------------------
        // Initial phase: prepare column and initialise first cell
        // ---------------------------------------------------------------------

        auto alignment_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();

        auto cell = *alignment_column_it;
        cell = this->track_cell(
                this->initialise_band_first_cell(cell.optimal_score(),
                                                 *++alignment_column_it,
                                                 this->m_scoring_scheme.score(sequence1_value,
                                                                              *std::ranges::begin(sequence2))),
                *cell_index_column_it);

        // ---------------------------------------------------------------------
        // Iteration phase: iterate over column and compute each cell
        // ---------------------------------------------------------------------

        for (auto && sequence2_value : sequence2 | views::drop(1))
        {
            auto cell = *alignment_column_it;
            cell = this->track_cell(
                this->compute_inner_cell(cell.optimal_score(),
                                         *++alignment_column_it,
                                         this->m_scoring_scheme.score(sequence1_value, sequence2_value)),
                *++cell_index_column_it);
        }

        // ---------------------------------------------------------------------
        // Final phase: track last cell
        // ---------------------------------------------------------------------

        this->track_last_row_cell(*alignment_column_it, *cell_index_column_it);
    }
};

} // namespace seqan3::detail
