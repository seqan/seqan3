// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_algorithm.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <memory>
#include <optional>
#include <type_traits>

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/configuration/align_config_result.hpp>
#include <seqan3/alignment/exception.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alignment/pairwise/detail/concept.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/alignment/matrix/detail/aligned_sequence_builder.hpp>

#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/view_to_simd.hpp>
#include <seqan3/core/type_traits/deferred_crtp_base.hpp>
#include <seqan3/core/type_traits/function.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/views/drop.hpp>
#include <seqan3/range/views/get.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief The alignment algorithm type to compute standard pairwise alignment using dynamic programming.
 * \implements std::invocable
 * \ingroup pairwise_alignment
 * \tparam config_t             The configuration type; must be of type seqan3::configuration.
 * \tparam algorithm_policies_t Template parameter pack with the policies to determine the execution of the algorithm;
 *                              must be wrapped as seqan3::detail::deferred_crtp_base.
 *
 * \details
 *
 * # Configuration
 *
 * The first template argument is the type of the alignment configuration which was used to configure the alignment
 * algorithm type. They must be the same, otherwise it is possible that the output is not coherent with the expected
 * result given the different configurations. The correct alignment is configured within the
 * seqan3::detail::alignment_configurator and returned as a std::function object, which can be passed around and
 * copied to create, for example multiple instances of the algorithm that can be executed in parallel.
 *
 * # Policies
 *
 * This type uses multiple inheritance to configure a specific alignment computation, e.g. global or local alignments,
 * using SIMD operations or scalar operations, computing the traceback or only the score etc.. These configurations
 * are inherited using so-called `alignment policies`. An alignment policy is a type that implements a specific
 * functionality through a common interface that is used by the alignment algorithm. These policies are also
 * the customisation points of the algorithm which will be used to implement a specific behaviour. You can read more
 * about the policies in \ref alignment_policy.
 *
 * Since some of the policies are augmented with traits to further refine the policy execution during the configuration,
 * it is necessary to defer the template instantiation of the policies, which are modelled as CRTP-base classes.
 * Therefore every added policy must be wrapped in a seqan3::detail::deferred_crtp_base class wrapper. This wrapper
 * then is expanded during the final template instantiation of the alignment algorithm type using the corresponding
 * template function seqan3::detail::invoke_deferred_crtp_base, which instantiates the CRTP policy types with the
 * correct alignment algorithm type.
 */
template <typename config_t, typename ...algorithm_policies_t>
class alignment_algorithm :
    public invoke_deferred_crtp_base<algorithm_policies_t, alignment_algorithm<config_t, algorithm_policies_t...>>...
{
private:

    //!\brief The alignment configuration traits type with auxiliary information extracted from the configuration type.
    using traits_t = alignment_configuration_traits<config_t>;
    //!\brief The type of an alignment column as defined by the respective matrix policy.
    using alignment_column_t = decltype(std::declval<alignment_algorithm>().current_alignment_column());
    //!\brief The iterator type over the alignment column.
    using alignment_column_iterator_t = std::ranges::iterator_t<alignment_column_t>;
    //!\brief The alignment result type.
    using alignment_result_t = typename traits_t::alignment_result_type;

    static_assert(!std::same_as<alignment_result_t, empty_type>, "Alignment result type was not configured.");

    //!\brief The type of the score debug matrix.
    using score_debug_matrix_t =
        std::conditional_t<traits_t::is_debug,
                           two_dimensional_matrix<std::optional<typename traits_t::original_score_type>,
                                                  std::allocator<std::optional<typename traits_t::original_score_type>>,
                                                  matrix_major_order::column>,
                           empty_type>;
    //!\brief The type of the trace debug matrix.
    using trace_debug_matrix_t =
        std::conditional_t<traits_t::is_debug,
                           two_dimensional_matrix<std::optional<trace_directions>,
                                                  std::allocator<std::optional<trace_directions>>,
                                                  matrix_major_order::column>,
                           empty_type>;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr alignment_algorithm()                                        = default; //!< Defaulted
    constexpr alignment_algorithm(alignment_algorithm const &)             = default; //!< Defaulted
    constexpr alignment_algorithm(alignment_algorithm &&)                  = default; //!< Defaulted
    constexpr alignment_algorithm & operator=(alignment_algorithm const &) = default; //!< Defaulted
    constexpr alignment_algorithm & operator=(alignment_algorithm &&)      = default; //!< Defaulted
    ~alignment_algorithm()                                                 = default; //!< Defaulted

    /*!\brief Constructs the algorithm with the passed configuration.
     * \param cfg The configuration to be passed to the algorithm.
     *
     * \details
     *
     * Maintains a copy of the configuration object on the heap using a std::shared_ptr. In addition, the alignment
     * state is initialised.
     */
    explicit constexpr alignment_algorithm(config_t const & cfg) : cfg_ptr{std::make_shared<config_t>(cfg)}
    {
        this->scoring_scheme = seqan3::get<align_cfg::scoring>(*cfg_ptr).value;
        this->initialise_alignment_state(*cfg_ptr);
    }
    //!\}

    /*!\name Invocation
     * \{
     */
    /*!\brief Computes the pairwise sequence alignment for the given range over indexed sequence pairs.
     * \tparam indexed_sequence_pairs_t The type of indexed_sequence_pairs; must model
     *                                  seqan3::detail::indexed_sequence_pair_range.
     * \tparam callback_t The type of the callback function that is called with the alignment result; must model
     *                    std::invocable accepting one argument of type seqan3::alignment_result.
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
        requires (!traits_t::is_vectorised) && std::invocable<callback_t, alignment_result_t>
    //!\endcond
    void operator()(indexed_sequence_pairs_t && indexed_sequence_pairs, callback_t && callback)
    {
        using std::get;

        for (auto && [sequence_pair, idx] : indexed_sequence_pairs)
            compute_single_pair(idx, get<0>(sequence_pair), get<1>(sequence_pair), callback);
    }

    //!\overload
    template <indexed_sequence_pair_range indexed_sequence_pairs_t, typename callback_t>
    //!\cond
        requires traits_t::is_vectorised && std::invocable<callback_t, alignment_result_t>
    //!\endcond
    void operator()(indexed_sequence_pairs_t && indexed_sequence_pairs, callback_t && callback)
    {
        assert(cfg_ptr != nullptr);

        static_assert(simd_concept<typename traits_t::score_type>, "Expected simd score type.");
        static_assert(simd_concept<typename traits_t::trace_type>, "Expected simd trace type.");

        // Extract the batch of sequences for the first and the second sequence.
        auto sequence1_range = indexed_sequence_pairs | views::get<0> | views::get<0>;
        auto sequence2_range = indexed_sequence_pairs | views::get<0> | views::get<1>;

        // Initialise the find_optimum policy in the simd case.
        this->initialise_find_optimum_policy(sequence1_range,
                                             sequence2_range,
                                             this->scoring_scheme.padding_match_score());

        // Convert batch of sequences to sequence of simd vectors.
        auto simd_sequences1 = convert_batch_of_sequences_to_simd_vector(sequence1_range);
        auto simd_sequences2 = convert_batch_of_sequences_to_simd_vector(sequence2_range);

        max_size_in_collection = std::pair{simd_sequences1.size(), simd_sequences2.size()};
        // Reset the alignment state's optimum between executions of the alignment algorithm.
        this->alignment_state.reset_optimum();

        compute_matrix(simd_sequences1, simd_sequences2);

        make_alignment_result(indexed_sequence_pairs, callback);
    }
    //!\}

private:
    /*!\brief Converts a batch of sequences to a sequence of simd vectors.
     * \tparam sequence_range_t The type of the range over sequences; must model std::ranges::forward_range.
     *
     * \param[in] sequences The batch of sequences to transform.
     *
     * \returns a sequence over simd vectors.
     *
     * \details
     *
     * Expects that the size of the batch is less or equal than the number of alignments that can be computed within one
     * simd vector. Applies an Array-of-Structures (AoS) to Structure-of-Arrays (SoA) transformation by storing one
     * column of the batch as a simd vector.
     */
    template <typename sequence_range_t>
    constexpr auto convert_batch_of_sequences_to_simd_vector(sequence_range_t & sequences)
    {
        assert(static_cast<size_t>(std::ranges::distance(sequences)) <= traits_t::alignments_per_vector);

        using simd_score_t = typename traits_t::score_type;

        std::vector<simd_score_t, aligned_allocator<simd_score_t, alignof(simd_score_t)>> simd_sequence{};

        for (auto && simd_vector_chunk : sequences | views::to_simd<simd_score_t>(traits_t::padding_symbol))
            for (auto && simd_vector : simd_vector_chunk)
                simd_sequence.push_back(std::move(simd_vector));

        return simd_sequence;
    }

    /*!\brief Computes the pairwise sequence alignment for a single pair of sequences.
     * \tparam sequence1_t The type of the first sequence; must model std::ranges::forward_range.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::forward_range.
     * \tparam callback_t The type of the callback function.
     *
     * \param[in] idx The index of the current processed sequence pair.
     * \param[in] sequence1 The first sequence (or packed sequences).
     * \param[in] sequence2 The second sequence (or packed sequences).
     * \param[in] callback The callback function to be invoked with the alignment result.
     *
     * \throws std::bad_alloc during allocation of the alignment matrices or
     *         seqan3::invalid_alignment_configuration if an invalid configuration for the given sequences is detected.
     *
     * \details
     *
     * Uses the standard dynamic programming algorithm to compute the pairwise sequence alignment.
     */
    template <std::ranges::forward_range sequence1_t,
              std::ranges::forward_range sequence2_t,
              typename callback_t>
    constexpr void compute_single_pair(size_t const idx,
                                       sequence1_t && sequence1,
                                       sequence2_t && sequence2,
                                       callback_t & callback)
    {
        assert(cfg_ptr != nullptr);

        if constexpr (traits_t::is_debug)
            initialise_debug_matrices(sequence1, sequence2);

        // Reset the alignment state's optimum between executions of the alignment algorithm.
        this->alignment_state.reset_optimum();

        if constexpr (traits_t::is_banded)
        {
            // Get the band and check if band configuration is valid.
            auto const & band = seqan3::get<align_cfg::band>(*cfg_ptr).value;
            check_valid_band_parameter(sequence1, sequence2, band);
            auto && [subsequence1, subsequence2] = this->slice_sequences(sequence1, sequence2, band);
            // It would be great to use this interface here instead
            compute_matrix(subsequence1, subsequence2, band);
            make_alignment_result(idx, subsequence1, subsequence2, callback);
        }
        else
        {
            compute_matrix(sequence1, sequence2);
            make_alignment_result(idx, sequence1, sequence2, callback);
        }
    }

    /*!\brief Checks if the band parameters are valid for the given sequences.
     * \tparam sequence1_t The type of the first sequence.
     * \tparam sequence2_t The type of the second sequence.
     *
     * \param[in] sequence1 The first sequence.
     * \param[in] sequence2 The second sequence.
     * \param[in] band The band to check.
     *
     * \throws seqan3::invalid_alignment_configuration if the band parameters would form an invalid alignment matrix.
     *
     * \details
     *
     * Checks if the given band intersects with the alignment matrix formed by the first and second sequence.
     * For example if the lower bound of the band is larger than the size of the first sequence the band would lie
     * outside of the alignment matrix and thus is invalid.
     */
    template <typename sequence1_t, typename sequence2_t>
    constexpr void check_valid_band_parameter(sequence1_t && sequence1,
                                              sequence2_t && sequence2,
                                              static_band const & band)
    {
        static_assert(config_t::template exists<align_cfg::band>(),
                      "The band configuration is required for the banded alignment algorithm.");

        using diff_type = std::iter_difference_t<std::ranges::iterator_t<sequence1_t>>;
        static_assert(std::is_signed_v<diff_type>,  "Only signed types can be used to test the band parameters.");

        if (static_cast<diff_type>(band.lower_bound) > std::ranges::distance(sequence1))
        {
            throw invalid_alignment_configuration
            {
                "Invalid band error: The lower bound excludes the whole alignment matrix."
            };
        }

        if (static_cast<diff_type>(band.upper_bound) < -std::ranges::distance(sequence2))
        {
            throw invalid_alignment_configuration
            {
                "Invalid band error: The upper bound excludes the whole alignment matrix."
            };
        }
    }

    /*!\brief Initialises the debug matrices for the given sequences.
     * \tparam sequence1_t The type of the first sequence.
     * \tparam sequence2_t The type of the second sequence.
     *
     * param[in] sequence1 The first sequence.
     * param[in] sequence2 The second sequence.
     *
     * \details
     *
     * Initialises the debug matrices if the alignment algorithm is running in debug mode. See seqan3::align_cfg::debug
     * for more information.
     */
    template <typename sequence1_t, typename sequence2_t>
    constexpr void initialise_debug_matrices(sequence1_t & sequence1, sequence2_t & sequence2)
    {
        size_t rows = std::ranges::distance(sequence2) + 1;
        size_t cols = std::ranges::distance(sequence1) + 1;

        score_debug_matrix = score_debug_matrix_t{number_rows{rows}, number_cols{cols}};
        trace_debug_matrix = trace_debug_matrix_t{number_rows{rows}, number_cols{cols}};
    }

    /*!\brief Compute the alignment by iterating over the alignment matrix in a column wise manner.
     * \tparam sequence1_t The type of the first sequence.
     * \tparam sequence2_t The type of the second sequence.
     *
     * \param[in] sequence1 The first sequence.
     * \param[in] sequence2 The second sequence.
     */
    template <typename sequence1_t, typename sequence2_t>
    void compute_matrix(sequence1_t & sequence1, sequence2_t & sequence2)
    //!\cond
        requires !traits_t::is_banded
    //!\endcond
    {
        // ----------------------------------------------------------------------------
        // Initialisation phase: allocate memory and initialise first column.
        // ----------------------------------------------------------------------------

        this->allocate_matrix(sequence1, sequence2);
        initialise_first_alignment_column(sequence2);

        // ----------------------------------------------------------------------------
        // Recursion phase: compute column-wise the alignment matrix.
        // ----------------------------------------------------------------------------

        for (auto const & seq1_value : sequence1)
        {
            compute_alignment_column<true>(seq1_value, sequence2);
            finalise_last_cell_in_column(true);
        }

        // ----------------------------------------------------------------------------
        // Wrap up phase: track score in last column and prepare the alignment result.
        // ----------------------------------------------------------------------------

        finalise_alignment();
    }

    //!\overload
    template <typename sequence1_t, typename sequence2_t>
    void compute_matrix(sequence1_t & sequence1, sequence2_t & sequence2, static_band const & band)
    //!\cond
        requires traits_t::is_banded
    //!\endcond
    {
        // ----------------------------------------------------------------------------
        // Initialisation phase: allocate memory and initialise first column.
        // ----------------------------------------------------------------------------

        // Allocate and initialise first column.
        this->allocate_matrix(sequence1, sequence2, band, this->alignment_state);
        size_t last_row_index = this->score_matrix.band_row_index;
        initialise_first_alignment_column(sequence2 | views::take(last_row_index));

        // ----------------------------------------------------------------------------
        // 1st recursion phase: iterate as long as the band intersects with the first row.
        // ----------------------------------------------------------------------------

        size_t sequence2_size = std::ranges::distance(sequence2);
        for (auto const & seq1_value : sequence1 | views::take(this->score_matrix.band_col_index))
        {
            compute_alignment_column<true>(seq1_value, sequence2 | views::take(++last_row_index));
            // Only if band reached last row of matrix the last cell might be tracked.
            finalise_last_cell_in_column(last_row_index >= sequence2_size);
        }

        // ----------------------------------------------------------------------------
        // 2nd recursion phase: iterate until the end of the matrix.
        // ----------------------------------------------------------------------------

        size_t first_row_index = 0;
        for (auto const & seq1_value : sequence1 | views::drop(this->score_matrix.band_col_index))
        {
            // In the second phase the band moves in every column one base down on the second sequence.
            compute_alignment_column<false>(seq1_value, sequence2 | views::slice(first_row_index++, ++last_row_index));
            // Only if band reached last row of matrix the last cell might be tracked.
            finalise_last_cell_in_column(last_row_index >= sequence2_size);
        }

        // ----------------------------------------------------------------------------
        // Wrap up phase: track score in last column and prepare the alignment result.
        // ----------------------------------------------------------------------------

        finalise_alignment();
    }

    /*!\brief Initialises the first column of the alignment matrix.
     * \tparam sequence2_t The type of the second sequence.
     *
     * \param[in] sequence2 The second sequence to initialise the first column for.
     *
     * \details
     *
     * Initialises the alignment matrix with special initialisation functions for the origin cell
     * and the cells in the first alignment column. The second sequence is required to get the size of the first
     * column which can vary between banded and unbanded alignments. The value of the second sequence is actually not
     * used during the initialisation.
     */
    template <typename sequence2_t>
    auto initialise_first_alignment_column(sequence2_t && sequence2)
    {
        // Get the initial column.
        alignment_column = this->current_alignment_column();
        assert(!alignment_column.empty()); // Must contain at least one element.

        // Initialise first cell.
        alignment_column_it = alignment_column.begin();
        this->init_origin_cell(*alignment_column_it, this->alignment_state);

        // Initialise the remaining cells of this column.
        for (auto it = std::ranges::begin(sequence2); it != std::ranges::end(sequence2); ++it)
            this->init_column_cell(*++alignment_column_it, this->alignment_state);

        // Finalise the last cell of the initial column.
        bool at_last_row = true;
        if constexpr (traits_t::is_banded) // If the band reaches until the last row of the matrix.
            at_last_row = static_cast<size_t>(this->score_matrix.band_row_index) == this->score_matrix.num_rows - 1;

        finalise_last_cell_in_column(at_last_row);
    }

    /*!\brief Computes a single alignment column.
     * \tparam initialise_first_cell An explicit bool template argument indicating whether the first cell of the current
     *                               alignment column needs to be initialised.
     * \tparam seq1_value_t The value type of the first sequence.
     * \tparam sequence2_t The type of the second sequence.
     *
     * \param[in] seq1_value The current value of the first sequence for this alignment column.
     * \param[in] sequence2 The current slice of the second sequence for this alignment column.
     *
     * \details
     *
     * Computes a single column within the alignment matrix. The first cell of the column is either initialised
     * according to the initialisation policy if `initialise_first_cell` is `true`, otherwise it assumes a banded
     * column within the matrix and computes the value accordingly.
     */
    template <bool initialise_first_cell, typename sequence1_value_t, typename sequence2_t>
    void compute_alignment_column(sequence1_value_t const & seq1_value, sequence2_t && sequence2)
    {
        this->next_alignment_column();  // move to next column and set alignment column iterator accordingly.
        alignment_column = this->current_alignment_column();
        alignment_column_it = alignment_column.begin();

        auto seq2_it = std::ranges::begin(sequence2);

        if constexpr (initialise_first_cell) // Initialise first cell if it intersects with the first row of the matrix.
        {
            this->init_row_cell(*alignment_column_it, this->alignment_state);
        }
        else // Compute first cell of banded column if it does not intersect with the first row of the matrix.
        {
            this->compute_first_band_cell(*alignment_column_it,
                                          this->alignment_state,
                                          this->scoring_scheme.score(seq1_value, *seq2_it));
            ++seq2_it;
        }

        for (; seq2_it != std::ranges::end(sequence2); ++seq2_it)
            this->compute_cell(*++alignment_column_it,
                               this->alignment_state,
                               this->scoring_scheme.score(seq1_value, *seq2_it));
    }

    /*!\brief Finalises the last cell of the current alignment column.
     * \param[in] at_last_row A bool indicating whether the column ends in the last row of the alignment matrix.
     *
     * \details
     *
     * After computing the alignment column the alignment column iterator possibly points to the last computed cell.
     * This value is used to check for a new optimum if searching in the last row was enabled. Otherwise it does
     * nothing. In case the alignment algorithm is run in debug mode the computed column is dumped in the debug
     * score and trace matrix.
     */
    constexpr void finalise_last_cell_in_column(bool const at_last_row) noexcept
    {
        if (at_last_row)
            this->check_score_of_last_row_cell(*alignment_column_it, this->alignment_state);

        if constexpr (traits_t::is_debug)
            dump_alignment_column();
    }

    //!\brief Checks the last cell, respectively column for the alignment optimum.
    constexpr void finalise_alignment() noexcept
    {
        // ----------------------------------------------------------------------------
        // Check for the optimum in last cell/column.
        // ----------------------------------------------------------------------------

        this->check_score_of_cells_in_last_column(alignment_column, this->alignment_state);
        this->check_score_of_last_cell(*alignment_column_it, this->alignment_state);
    }

    /*!\brief Creates a new alignment result from the current alignment optimum and for the given pair of sequences.
     * \tparam callback_t The type of the callback function.
     * \tparam index_t The type of the index.
     * \tparam sequence1_t The type of the first sequence.
     * \tparam sequence2_t The type of the second sequence.
     *
     * \param[in] idx The internal index used for this pair of sequences.
     * \param[in] sequence1 The first range to get the alignment for if requested.
     * \param[in] sequence2 The second range to get the alignment for if requested.
     * \param[in] callback The callback function to be invoked with the alignment result.
     *
     * \details
     *
     * Fills the alignment result with the requested values. Depending on the selected configuration the following
     * is extracted and/or computed:
     *
     * 1. The alignment score.
     * 2. The end positions of the aligned range for the first and second sequence.
     * 3. The begin positions of the aligned range for the first and second sequence.
     * 4. The alignment between both sequences in the respective aligned region.
     *
     * If the alignment is run in debug mode (see seqan3::align_cfg::debug) the debug score and optionally trace matrix
     * are stored in the alignment result as well.
     *
     * Finally, the callback is invoked with the computed alignment result.
     */
    template <typename index_t, typename sequence1_t, typename sequence2_t, typename callback_t>
    //!\cond
        requires !traits_t::is_vectorised
    //!\endcond
    constexpr void make_alignment_result(index_t const idx,
                                         sequence1_t & sequence1,
                                         sequence2_t & sequence2,
                                         callback_t & callback)
    {
        using result_value_t = typename alignment_result_value_type_accessor<alignment_result_t>::type;

        // ----------------------------------------------------------------------------
        // Build the alignment result
        // ----------------------------------------------------------------------------

        static_assert(config_t::template exists<align_cfg::result>(),
                      "The configuration must contain an align_cfg::result element.");

        result_value_t res{};

        res.id = idx;

        // Choose what needs to be computed.
        if constexpr (traits_t::compute_score)
            res.score = this->alignment_state.optimum.score;

        if constexpr (traits_t::compute_back_coordinate)
        {
            res.back_coordinate = alignment_coordinate{column_index_type{this->alignment_state.optimum.column_index},
                                                       row_index_type{this->alignment_state.optimum.row_index}};
            // At some point this needs to be refactored so that it is not necessary to adapt the coordinate.
            if constexpr (traits_t::is_banded)
                res.back_coordinate.second += res.back_coordinate.first - this->trace_matrix.band_col_index;
        }

        if constexpr (traits_t::compute_front_coordinate)
        {
            // Get a aligned sequence builder for banded or un-banded case.
            aligned_sequence_builder builder{sequence1, sequence2};
            auto optimum_coordinate = alignment_coordinate{column_index_type{this->alignment_state.optimum.column_index},
                                                           row_index_type{this->alignment_state.optimum.row_index}};
            auto trace_res = builder(this->trace_matrix.trace_path(optimum_coordinate));
            res.front_coordinate.first = trace_res.first_sequence_slice_positions.first;
            res.front_coordinate.second = trace_res.second_sequence_slice_positions.first;

            if constexpr (traits_t::compute_sequence_alignment)
                res.alignment = std::move(trace_res.alignment);
        }

        // Store the matrices in debug mode.
        if constexpr (traits_t::is_debug)
        {
            res.score_debug_matrix = std::move(score_debug_matrix);
            if constexpr (traits_t::result_type_rank == 3) // compute alignment
                res.trace_debug_matrix = std::move(trace_debug_matrix);
        }

        callback(std::move(res));
    }

    /*!\brief Creates a new alignment result from the current alignment optimum and for the given indexed sequence
     *        range.
     * \tparam callback_t The type of the callback function.
     * \tparam indexed_sequence_pair_range_t The type of the indexed sequence pair range.
     *
     * \param[in] index_sequence_pairs The range over indexed sequence pairs.
     * \param[in] callback The callback function to be invoked with the alignment result.
     *
     * \details
     *
     * This function is called for the vectorised algorithm. In this case the alignment state stores the results for
     * the entire chunk of sequence pairs processed within this alignment computation. Accordingly, the chunk of
     * sequence pairs is processed iteratively and the alignment results are added to the returned vector.
     * Depending on the selected configuration the following is extracted and/or computed:
     *
     * 1. The alignment score.
     * 2. The end positions of the aligned range for the first and second sequence.
     * 3. The begin positions of the aligned range for the first and second sequence.
     * 4. The alignment between both sequences in the respective aligned region.
     *
     * If the alignment is run in debug mode (see seqan3::align_cfg::debug) the debug score and optionally trace matrix
     * are stored in the alignment result as well.
     *
     * Finally, the callback is invoked with each computed alignment result iteratively.
     */
    template <typename indexed_sequence_pair_range_t, typename callback_t>
    //!\cond
        requires traits_t::is_vectorised
    //!\endcond
    constexpr auto make_alignment_result(indexed_sequence_pair_range_t && index_sequence_pairs,
                                         callback_t & callback)
    {
        using result_value_t = typename alignment_result_value_type_accessor<alignment_result_t>::type;

        size_t simd_index = 0;
        for (auto && [sequence_pairs, alignment_index] : index_sequence_pairs)
        {
            (void) sequence_pairs;
            result_value_t res{};
            res.id = alignment_index;

            if constexpr (traits_t::compute_score)
                res.score = this->alignment_state.optimum.score[simd_index];  // Just take this

            if constexpr (traits_t::compute_back_coordinate)
            {
                res.back_coordinate.first = this->alignment_state.optimum.column_index[simd_index] ;
                res.back_coordinate.second = this->alignment_state.optimum.row_index[simd_index];
            }

            callback(std::move(res));
            ++simd_index;
        }
    }

    /*!\brief Dumps the current alignment matrix in the debug score matrix and if requested debug trace matrix.
     *
     * \details
     *
     * Copies the current score and if configured the current trace column into the local debug score and trace matrix.
     * In the banded case the full matrix will be allocated with std::optional values and only the respective parts
     * of the matrix that are present in the band are filled.
     */
    void dump_alignment_column()
    {
        using std::get;

        auto column = this->current_alignment_column();

        auto coord = get<1>(column.front()).coordinate;
        if constexpr (traits_t::is_banded)
            coord.second += coord.first - this->score_matrix.band_col_index;

        matrix_offset offset{row_index_type{static_cast<std::ptrdiff_t>(coord.second)},
                             column_index_type{static_cast<std::ptrdiff_t>(coord.first)}};

        std::ranges::copy(column | std::views::transform([] (auto const & tpl)
        {
            using std::get;
            return get<0>(tpl).current;
        }), score_debug_matrix.begin() + offset);

        // if traceback is enabled.
        if constexpr (config_t::template exists<align_cfg::result<with_alignment_type>>())
        {
            auto trace_matrix_it = trace_debug_matrix.begin() + offset;
            std::ranges::copy(column | std::views::transform([] (auto const & tpl)
            {
                using std::get;
                return get<1>(tpl).current;
            }), trace_debug_matrix.begin() + offset);
        }
    }

    //!\brief The alignment configuration stored on the heap.
    std::shared_ptr<config_t> cfg_ptr{};
    //!\brief Stores the currently processed alignment column.
    alignment_column_t alignment_column{};
    //!\brief Stores the state of the currently processed alignment column.
    alignment_column_iterator_t alignment_column_it{};
    //!\brief The debug matrix for the scores.
    score_debug_matrix_t score_debug_matrix{};
    //!\brief The debug matrix for the traces.
    trace_debug_matrix_t trace_debug_matrix{};
    //!\brief The maximal size within the first and the second sequence collection.
    std::pair<size_t, size_t> max_size_in_collection{};
};

} // namespace seqan3::detail
