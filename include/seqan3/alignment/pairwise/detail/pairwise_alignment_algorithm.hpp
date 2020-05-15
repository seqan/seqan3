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

#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/matrix/detail/coordinate_matrix.hpp>
#include <seqan3/alignment/matrix/detail/score_matrix_single_column.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/detail/type_inspection.hpp>
#include <seqan3/range/views/drop.hpp>
#include <seqan3/range/views/zip.hpp>

namespace seqan3::detail
{

/*!\brief The alignment algorithm type to compute standard pairwise alignment using dynamic programming.
 * \implements std::invocable
 * \ingroup pairwise_alignment
 *
 * \tparam alignment_configuration_t The configuration type; must be of type seqan3::configuration.
 * \tparam policies_t Variadic template argument for the different policies of this alignment algorithm.
 *
 * \details
 *
 * ### Configuration
 *
 * The first template argument is the type of the alignment configuration. The alignment configuration was used to
 * configure the `alignment algorithm type` within the seqan3::detail::alignment_configurator.
 * The algorithm computes a column based dynamic programming matrix given two sequences.
 * After the computation a user defined callback function is invoked with the computed seqan3::alignment_result.
 */
template <typename alignment_configuration_t, typename ...policies_t>
//!\cond
    requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
//!\endcond
class pairwise_alignment_algorithm : protected policies_t...
{
private:
    //!\brief The alignment configuration traits type with auxiliary information extracted from the configuration type.
    using traits_type = alignment_configuration_traits<alignment_configuration_t>;
    //!\brief The type of the scoring scheme.
    using scoring_scheme_type =  typename traits_type::scoring_scheme_type;
    //!\brief The configured alignment result type.
    using alignment_result_type = typename traits_type::alignment_result_type;

    static_assert(!std::same_as<alignment_result_type, empty_type>, "Alignment result type was not configured.");

    //!\brief The configured scoring scheme.
    scoring_scheme_type m_scoring_scheme{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    pairwise_alignment_algorithm() = default; //!< Defaulted.
    pairwise_alignment_algorithm(pairwise_alignment_algorithm const &) = default; //!< Defaulted.
    pairwise_alignment_algorithm(pairwise_alignment_algorithm &&) = default; //!< Defaulted.
    pairwise_alignment_algorithm & operator=(pairwise_alignment_algorithm const &) = default; //!< Defaulted.
    pairwise_alignment_algorithm & operator=(pairwise_alignment_algorithm &&) = default; //!< Defaulted.
    ~pairwise_alignment_algorithm() = default; //!< Defaulted.

    /*!\brief Constructs and initialises the algorithm using the alignment configuration.
     * \param config The configuration passed into the algorithm.
     *
     * \details
     *
     * Initialises the algorithm given the user settings from the alignment configuration object.
     */
    pairwise_alignment_algorithm(alignment_configuration_t const & config) : policies_t(config)...
    {
        m_scoring_scheme = seqan3::get<align_cfg::scoring>(config).value;
    }
    //!\}

    /*!\name Invocation
     * \{
     */
    /*!\brief Computes the pairwise sequence alignment for the given range over indexed sequence pairs.
     * \tparam indexed_sequence_pairs_t The type of indexed_sequence_pairs; must model
     *                                  seqan3::detail::indexed_sequence_pair_range.
     * \tparam callback_t The type of the callback function that is called with the alignment result; must model
    *                    std::invocable with seqan3::alignment_result as argument.
     *
     * \param[in] indexed_sequence_pairs A range over indexed sequence pairs to be aligned.
     * \param[in] callback The callback function to be invoked with each computed alignment result.
     *
     * \throws std::bad_alloc during allocation of the alignment matrices or
     *         seqan3::invalid_alignment_configuration if an invalid configuration for the given sequences is detected.
     *
     * \details
     *
     * Uses the standard dynamic programming algorithm to compute the pairwise sequence alignment for each
     * sequence pair. The space and runtime complexities depend on the selected configurations (see below).
     * For every computed alignment the given callback is invoked with the respective alignment result.
     *
     * ### Exception
     *
     * Strong exception guarantee. Might throw std::bad_alloc or seqan3::invalid_alignment_configuration.
     *
     * ### Thread-safety
     *
     * Calls to this functions in a concurrent environment are not thread safe. Instead use a copy of the alignment
     * algorithm type.
     *
     * ### Complexity
     *
     * The following table lists the runtime and space complexities for the banded and unbanded algorithm dependent
     * on the configured seqan3::align_cfg::result per sequence pair.
     * Let `n` be the length of the first sequence, `m` be the length of the second sequence and `k` be the size of
     * the band.
     *
     * |                        | unbanded         | banded            |
     * |:----------------------:|:----------------:|:-----------------:|
     * |runtime                 |\f$ O(n*m) \f$    |\f$ O(n*k) \f$     |
     * |space (score only)      |\f$ O(m) \f$      |\f$ O(k) \f$       |
     * |space (end positions)   |\f$ O(m) \f$      |\f$ O(k) \f$       |
     * |space (begin positions) |\f$ O(n*m) \f$    |\f$ O(n*k) \f$     |
     * |space (alignment)       |\f$ O(n*m) \f$    |\f$ O(n*k) \f$     |
     */
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
            res.score = compute_matrix(get<0>(sequence_pair), get<1>(sequence_pair));
            callback(alignment_result_type{res});
        }
    }

protected:
    /*!\brief Compute the actual alignment.
     * \tparam sequence1_t The type of the first sequence; must model std::ranges::forward_range.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::forward_range.
     *
     * \param[in] sequence1 The first sequence to compute the alignment for.
     * \param[in] sequence2 The second sequence to compute the alignment for.
     */
    template <std::ranges::forward_range sequence1_t, std::ranges::forward_range sequence2_t>
    int32_t compute_matrix(sequence1_t && sequence1, sequence2_t && sequence2)
    {
        // ---------------------------------------------------------------------
        // Initialisation phase: allocate memory and initialise first column.
        // ---------------------------------------------------------------------

        this->reset_optimum(); // Reset the tracker for the new alignment computation.

        thread_local score_matrix_single_column<int32_t> local_score_matrix{};
        coordinate_matrix local_index_matrix{};

        size_t number_of_columns = std::ranges::distance(sequence1) + 1;
        size_t number_of_rows = std::ranges::distance(sequence2) + 1;

        local_score_matrix.resize(column_index_type{number_of_columns}, row_index_type{number_of_rows});
        local_index_matrix.resize(column_index_type{number_of_columns}, row_index_type{number_of_rows});

        auto alignment_matrix_it = local_score_matrix.begin();
        auto indexed_matrix_it = local_index_matrix.begin();

        initialise_column(*alignment_matrix_it, *indexed_matrix_it, sequence2);

        // ---------------------------------------------------------------------
        // Iteration phase: compute column-wise the alignment matrix.
        // ---------------------------------------------------------------------

        for (auto sequence1_value : sequence1)
            compute_column(*++alignment_matrix_it, *++indexed_matrix_it, sequence1_value, sequence2);

        // ---------------------------------------------------------------------
        // Final phase: track score of last column
        // ---------------------------------------------------------------------

        auto && alignment_column = *alignment_matrix_it;
        auto && cell_index_column = *indexed_matrix_it;

        auto alignment_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();

        this->track_last_column_cell(*alignment_column_it, *cell_index_column_it);

        for ([[maybe_unused]] auto && unused : sequence2)
            this->track_last_column_cell(*++alignment_column_it, *++cell_index_column_it);

        this->track_final_cell(*alignment_column_it, *cell_index_column_it);

        return this->tracked_optimum();
    }

    /*!\brief Initialise the first column of the alignment matrix.
     * \tparam alignment_column_t The type of the alignment column; must model std::ranges::input_range.
     * \tparam cell_index_column_t The type of the indexed column; must model std::ranges::input_range.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::input_range.
     *
     * \param[in] alignment_column The current alignment matrix column to compute.
     * \param[in] cell_index_column The current index matrix column to get the respective cell indices.
     * \param[in] sequence2 The second sequence used to determine the size of the column.
     *
     * \details
     *
     * The first column of the alignment matrix does not require any character comparisons of the sequences that
     * shall be aligned. The second sequence is thus only needed to determine the size of the column.
     * The computation of the column is split into three phases: the initialisation phase, the iteration phase, and
     * the final phase. In the initialisation phase the first cell of the column is computed and in the iteration
     * phase all remaining cells are computed. In the final phase the last cell is possibly evaluated for a new
     * alignment optimum.
     */
    template <std::ranges::input_range alignment_column_t,
              std::ranges::input_range cell_index_column_t,
              std::ranges::input_range sequence2_t>
    void initialise_column(alignment_column_t && alignment_column,
                           cell_index_column_t && cell_index_column,
                           sequence2_t && sequence2)
    {
        // ---------------------------------------------------------------------
        // Initial phase: prepare column and initialise first cell
        // ---------------------------------------------------------------------

        auto first_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();
        *first_column_it = this->track_cell(this->initialise_origin_cell(), *cell_index_column_it);

        // ---------------------------------------------------------------------
        // Iteration phase: iterate over column and compute each cell
        // ---------------------------------------------------------------------

        for ([[maybe_unused]] auto && unused : sequence2)
        {
            ++first_column_it;
            *first_column_it = this->track_cell(this->initialise_first_column_cell(*first_column_it),
                                                *++cell_index_column_it);
        }

        // ---------------------------------------------------------------------
        // Final phase: track last cell of initial column
        // ---------------------------------------------------------------------

        this->track_last_row_cell(*first_column_it, *cell_index_column_it);
    }

    /*!\brief Initialise any column of the alignment matrix except the first one.
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
     */
    template <std::ranges::input_range alignment_column_t,
              std::ranges::input_range cell_index_column_t,
              semialphabet sequence1_value_t,
              std::ranges::input_range sequence2_t>
    void compute_column(alignment_column_t && alignment_column,
                        cell_index_column_t && cell_index_column,
                        sequence1_value_t const & sequence1_value,
                        sequence2_t && sequence2)
    {
        using score_type = typename traits_type::score_type;

        // ---------------------------------------------------------------------
        // Initial phase: prepare column and initialise first cell
        // ---------------------------------------------------------------------

        auto alignment_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();

        auto cell = *alignment_column_it;
        score_type diagonal = cell.optimal_score();
        *alignment_column_it = this->track_cell(this->initialise_first_row_cell(cell), *cell_index_column_it);

        // ---------------------------------------------------------------------
        // Iteration phase: iterate over column and compute each cell
        // ---------------------------------------------------------------------

        for (auto && sequence2_value : sequence2)
        {
            auto cell = *++alignment_column_it;
            score_type next_diagonal = cell.optimal_score();
            *alignment_column_it = this->track_cell(
                this->compute_inner_cell(diagonal, cell, m_scoring_scheme.score(sequence1_value, sequence2_value)),
                *++cell_index_column_it);
            diagonal = next_diagonal;
        }

        // ---------------------------------------------------------------------
        // Final phase: track last cell
        // ---------------------------------------------------------------------

        this->track_last_row_cell(*alignment_column_it, *cell_index_column_it);
    }
};

} // namespace seqan3::detail
