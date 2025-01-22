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

#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/utility/container/aligned_allocator.hpp>
#include <seqan3/utility/detail/type_name_as_string.hpp>
#include <seqan3/utility/simd/views/to_simd.hpp>
#include <seqan3/utility/views/elements.hpp>

namespace seqan3::detail
{

/*!\brief The alignment algorithm type to compute standard pairwise alignment using dynamic programming.
 * \implements std::invocable
 * \ingroup alignment_pairwise
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
template <typename alignment_configuration_t, typename... policies_t>
    requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
class pairwise_alignment_algorithm : protected policies_t...
{
protected:
    //!\brief The alignment configuration traits type with auxiliary information extracted from the configuration type.
    using traits_type = alignment_configuration_traits<alignment_configuration_t>;
    //!\brief The configured score type.
    using score_type = typename traits_type::score_type;
    //!\brief The configured alignment result type.
    using alignment_result_type = typename traits_type::alignment_result_type;

    static_assert(!std::same_as<alignment_result_type, empty_type>, "Alignment result type was not configured.");

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    pairwise_alignment_algorithm() = default;                                                 //!< Defaulted.
    pairwise_alignment_algorithm(pairwise_alignment_algorithm const &) = default;             //!< Defaulted.
    pairwise_alignment_algorithm(pairwise_alignment_algorithm &&) = default;                  //!< Defaulted.
    pairwise_alignment_algorithm & operator=(pairwise_alignment_algorithm const &) = default; //!< Defaulted.
    pairwise_alignment_algorithm & operator=(pairwise_alignment_algorithm &&) = default;      //!< Defaulted.
    ~pairwise_alignment_algorithm() = default;                                                //!< Defaulted.

    /*!\brief Constructs and initialises the algorithm using the alignment configuration.
     * \param config The configuration passed into the algorithm.
     *
     * \details
     *
     * Initialises the base policies of the alignment algorithm.
     */
    pairwise_alignment_algorithm(alignment_configuration_t const & config) : policies_t(config)...
    {}
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
     * on the given \ref seqan3_align_cfg_output_configurations "seqan3::align_cfg::output_*" per sequence pair.
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
        requires std::invocable<callback_t, alignment_result_type>
    void operator()(indexed_sequence_pairs_t && indexed_sequence_pairs, callback_t && callback)
    {
        using std::get;

        for (auto && [sequence_pair, idx] : indexed_sequence_pairs)
        {
            size_t const sequence1_size = std::ranges::distance(get<0>(sequence_pair));
            size_t const sequence2_size = std::ranges::distance(get<1>(sequence_pair));

            auto && [alignment_matrix, index_matrix] = this->acquire_matrices(sequence1_size, sequence2_size);

            compute_matrix(get<0>(sequence_pair), get<1>(sequence_pair), alignment_matrix, index_matrix);
            this->make_result_and_invoke(std::forward<decltype(sequence_pair)>(sequence_pair),
                                         std::move(idx),
                                         this->optimal_score,
                                         this->optimal_coordinate,
                                         alignment_matrix,
                                         callback);
        }
    }
    //!\}

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

        convert_batch_of_sequences_to_simd_vector(simd_seq1_collection,
                                                  seq1_collection,
                                                  this->scoring_scheme.padding_symbol);
        convert_batch_of_sequences_to_simd_vector(simd_seq2_collection,
                                                  seq2_collection,
                                                  this->scoring_scheme.padding_symbol);

        size_t const sequence1_size = std::ranges::distance(simd_seq1_collection);
        size_t const sequence2_size = std::ranges::distance(simd_seq2_collection);

        auto && [alignment_matrix, index_matrix] = this->acquire_matrices(sequence1_size, sequence2_size);

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

protected:
    /*!\brief Converts a batch of sequences to a sequence of simd vectors.
     * \tparam simd_sequence_t The type of the simd sequence; must model std::ranges::output_range for the `score_type`.
     * \tparam sequence_collection_t The type of the collection containing the sequences; must model
     *                               std::ranges::forward_range.
     * \tparam padding_symbol_t The type of the padding symbol.
     *
     * \param[out] simd_sequence The transformed simd sequence.
     * \param[in] sequences The batch of sequences to transform.
     * \param[in] padding_symbol The symbol that should be appended during the transformation/
     *
     * \details
     *
     * Expects that the size of the collection is less or equal than the number of alignments that can be computed
     * within one simd vector (simd_traits\<score_type\>\::length).
     * Applies an Array-of-Structures (AoS) to Structure-of-Arrays (SoA) transformation by storing one
     * column of the collection as a simd vector. The resulting simd sequence has the size of the longest sequence in
     * the collection. For all sequences with a smaller size the padding symbol will be appended during the simd
     * transformation to fill up the remaining size difference.
     */
    template <typename simd_sequence_t, std::ranges::forward_range sequence_collection_t, arithmetic padding_symbol_t>
        requires std::ranges::output_range<simd_sequence_t, score_type>
    void convert_batch_of_sequences_to_simd_vector(simd_sequence_t & simd_sequence,
                                                   sequence_collection_t & sequences,
                                                   padding_symbol_t const & padding_symbol)
    {
        assert(static_cast<size_t>(std::ranges::distance(sequences)) <= traits_type::alignments_per_vector);

        simd_sequence.clear();
        for (auto && simd_vector_chunk : sequences | views::to_simd<score_type>(padding_symbol))
            std::ranges::move(simd_vector_chunk, std::back_inserter(simd_sequence));
    }

    /*!\brief Compute the actual alignment.
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

        initialise_column(*alignment_matrix_it, *indexed_matrix_it, sequence2);

        // ---------------------------------------------------------------------
        // Iteration phase: compute column-wise the alignment matrix.
        // ---------------------------------------------------------------------

        for (auto alphabet1 : sequence1)
            compute_column(*++alignment_matrix_it,
                           *++indexed_matrix_it,
                           this->scoring_scheme_profile_column(alphabet1),
                           sequence2);

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

        for ([[maybe_unused]] auto const & unused : sequence2)
        {
            ++first_column_it;
            *first_column_it =
                this->track_cell(this->initialise_first_column_cell(*first_column_it), *++cell_index_column_it);
        }

        // ---------------------------------------------------------------------
        // Final phase: track last cell of initial column
        // ---------------------------------------------------------------------

        this->track_last_row_cell(*first_column_it, *cell_index_column_it);
    }

    /*!\brief Initialise any column of the alignment matrix except the first one.
     * \tparam alignment_column_t The type of the alignment column; must model std::ranges::input_range.
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
     */
    template <std::ranges::input_range alignment_column_t,
              std::ranges::input_range cell_index_column_t,
              typename alphabet1_t,
              std::ranges::input_range sequence2_t>
        requires semialphabet<alphabet1_t> || simd_concept<alphabet1_t>
    void compute_column(alignment_column_t && alignment_column,
                        cell_index_column_t && cell_index_column,
                        alphabet1_t const & alphabet1,
                        sequence2_t && sequence2)
    {
        using score_type = typename traits_type::score_type;

        // ---------------------------------------------------------------------
        // Initial phase: prepare column and initialise first cell
        // ---------------------------------------------------------------------

        auto alignment_column_it = alignment_column.begin();
        auto cell_index_column_it = cell_index_column.begin();

        auto cell = *alignment_column_it;
        score_type diagonal = cell.best_score();
        *alignment_column_it = this->track_cell(this->initialise_first_row_cell(cell), *cell_index_column_it);

        // ---------------------------------------------------------------------
        // Iteration phase: iterate over column and compute each cell
        // ---------------------------------------------------------------------

        for (auto const & alphabet2 : sequence2)
        {
            auto cell = *++alignment_column_it;
            score_type next_diagonal = cell.best_score();
            *alignment_column_it = this->track_cell(
                this->compute_inner_cell(diagonal, cell, this->scoring_scheme.score(alphabet1, alphabet2)),
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
