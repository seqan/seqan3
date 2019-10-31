// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
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
#include <seqan3/alignment/matrix/detail/aligned_sequence_builder.hpp>

#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/type_traits/deferred_crtp_base.hpp>
#include <seqan3/range/views/drop.hpp>
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
    //!\brief Check if the alignment is banded.
    static constexpr bool is_banded = config_t::template exists<align_cfg::band>();
    //!\brief Check if debug mode is enabled.
    static constexpr bool is_debug_mode = config_t::template exists<detail::algorithm_debugging>();

    //!\brief The type of the stored scoring scheme.
    using score_scheme_t = decltype(seqan3::get<align_cfg::scoring>(std::declval<config_t>()).value);
    //!\brief The type of a column inside of the score matrix.
    using score_matrix_column_t =
        std::ranges::iter_value_t<decltype(std::declval<alignment_algorithm>().score_matrix_iter)>;
    //!\brief The type of a column inside of the trace matrix.
    using trace_matrix_column_t =
        std::ranges::iter_value_t<decltype(std::declval<alignment_algorithm>().trace_matrix_iter)>;
    //!\brief The type of an alignment column as defined by the respective matrix policy.
    using alignment_column_t = decltype(std::declval<alignment_algorithm>().current_alignment_column());
    //!\brief The iterator type over the alignment column.
    using alignment_column_iterator_t = std::ranges::iterator_t<alignment_column_t>;

    //!\brief The type of the score debug matrix.
    using score_debug_matrix_t =
        std::conditional_t<is_debug_mode,
                           two_dimensional_matrix<std::optional<int32_t>,
                                                  std::allocator<std::optional<int32_t>>,
                                                  matrix_major_order::column>,
                           empty_type>;
    //!\brief The type of the trace debug matrix.
    using trace_debug_matrix_t =
        std::conditional_t<is_debug_mode,
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
        score_scheme = seqan3::get<align_cfg::scoring>(*cfg_ptr).value;
        this->initialise_alignment_state(*cfg_ptr);
    }
    //!\}

    /*!\name Invocation
     * \{
     */
    /*!\brief Computes the pairwise sequence alignment for the given pair of sequences.
     * \tparam sequence1_t The type of the first sequence; must model std::ranges::forward_range.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::forward_range.
     *
     * \param[in] idx The index of the current processed sequence pair.
     * \param[in] sequence1 The first sequence (or packed sequences).
     * \param[in] sequence2 The second sequence (or packed sequences).
     *
     * \returns A seqan3::alignment_result with the requested alignment outcomes.
     *
     * \throws std::bad_alloc during allocation of the alignment matrices or
     *         seqan3::invalid_alignment_configuration if an invalid configuration for the given sequences is detected.
     *
     * \details
     *
     * Uses the standard dynamic programming algorithm to compute the pairwise sequence alignment. The space and
     * runtime complexities depend on the selected configurations (see below).
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
     * on the configured seqan3::align_cfg::result.
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
    template <std::ranges::forward_range sequence1_t, std::ranges::forward_range sequence2_t>
    auto operator()(size_t const idx, sequence1_t && sequence1, sequence2_t && sequence2)
    {
        assert(cfg_ptr != nullptr);

        if constexpr (is_debug_mode)
            initialise_debug_matrices(sequence1, sequence2);

        // Reset the alignment state's optimum between executions of the alignment algorithm.
        this->alignment_state.reset_optimum();

        return compute_matrix(idx, sequence1, sequence2);
    }
    //!\}

private:
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
     * \tparam index_t The type of the index.
     * \tparam sequence1_t The type of the first sequence.
     * \tparam sequence2_t The type of the second sequence.
     *
     * \param[in] idx The corresponding index for the given sequence pair.
     * \param[in] sequence1 The first sequence.
     * \param[in] sequence2 The second sequence.
     */
    template <typename index_t, typename sequence1_t, typename sequence2_t>
    auto compute_matrix(index_t const idx, sequence1_t & sequence1, sequence2_t & sequence2)
    //!\cond
        requires !is_banded
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

        return finalise_alignment(idx, sequence1, sequence2);
    }

    //!\overload
    template <typename index_t, typename sequence1_t, typename sequence2_t>
    auto compute_matrix(index_t const idx, sequence1_t & sequence1, sequence2_t & sequence2)
    //!\cond
        requires is_banded
    //!\endcond
    {
        // ----------------------------------------------------------------------------
        // Initialisation phase: allocate memory and initialise first column.
        // ----------------------------------------------------------------------------

        // Get the band and check if band configuration is valid.
        auto const & band = seqan3::get<align_cfg::band>(*cfg_ptr).value;
        check_valid_band_parameter(sequence1, sequence2, band);

        // Slice sequences such that band starts in origin and ends in sink.
        auto && [seq1_slice, seq2_slice] = this->slice_sequences(sequence1, sequence2, band);

        // Allocate and initialise first column.
        this->allocate_matrix(seq1_slice, seq2_slice, band, this->alignment_state);
        size_t last_row_index = this->score_matrix.band_row_index;
        initialise_first_alignment_column(seq2_slice | views::take(last_row_index));

        // ----------------------------------------------------------------------------
        // 1st recursion phase: iterate as long as the band intersects with the first row.
        // ----------------------------------------------------------------------------

        size_t seq2_slice_size = std::ranges::distance(seq2_slice);
        for (auto const & seq1_value : seq1_slice | views::take(this->score_matrix.band_col_index))
        {
            compute_alignment_column<true>(seq1_value, seq2_slice | views::take(++last_row_index));
            // Only if band reached last row of matrix the last cell might be tracked.
            finalise_last_cell_in_column(last_row_index >= seq2_slice_size);
        }

        // ----------------------------------------------------------------------------
        // 2nd recursion phase: iterate until the end of the matrix.
        // ----------------------------------------------------------------------------

        size_t first_row_index = 0;
        for (auto const & seq1_value : seq1_slice | views::drop(this->score_matrix.band_col_index))
        {
            // In the second phase the band moves in every column one base down on the second sequence.
            compute_alignment_column<false>(seq1_value, seq2_slice | views::slice(first_row_index++, ++last_row_index));
            // Only if band reached last row of matrix the last cell might be tracked.
            finalise_last_cell_in_column(last_row_index >= seq2_slice_size);
        }

        // ----------------------------------------------------------------------------
        // Wrap up phase: track score in last column and prepare the alignment result.
        // ----------------------------------------------------------------------------

        return finalise_alignment(idx, seq1_slice, seq2_slice);
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
        if constexpr (is_banded) // If the band reaches until the last row of the matrix.
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
                                          score_scheme.score(seq1_value, *seq2_it));
            ++seq2_it;
        }

        for (; seq2_it != std::ranges::end(sequence2); ++seq2_it)
            this->compute_cell(*++alignment_column_it, this->alignment_state, score_scheme.score(seq1_value, *seq2_it));
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

        if constexpr (is_debug_mode)
            dump_alignment_column();
    }

    /*!\brief Creates a new alignment result from the current alignment optimum and for the given pair of sequences.
     * \tparam index_t The type of the index.
     * \tparam sequence1_t The type of the first sequence.
     * \tparam sequence2_t The type of the second sequence.
     *
     * \param[in] idx The internal index used for this pair of sequences.
     * \param[in] sequence1 The first range to get the alignment for if requested.
     * \param[in] sequence2 The second range to get the alignment for if requested.
     *
     * \returns A seqan3::alignment_result with the requested alignment outcomes.
     *
     * \details
     *
     * At first the last column/cell of the alignment matrix is checked for a new alignment optimum.
     * Then the alignment result is prepared. Depending on the selected configuration the following is extracted and/or
     * computed:
     *
     * 1. The alignment score.
     * 2. The end positions of the aligned range for the first and second sequence.
     * 3. The begin positions of the aligned range for the first and second sequence.
     * 4. The alignment between both sequences in the respective aligned region.
     *
     * If the alignment is run in debug mode (see seqan3::align_cfg::debug) the debug score and optionally trace matrix
     * are stored in the alignment result as well.
     */
    template <typename index_t, typename sequence1_t, typename sequence2_t>
    constexpr auto finalise_alignment(index_t const & idx,
                                      sequence1_t && sequence1,
                                      sequence2_t && sequence2)
    {
        // ----------------------------------------------------------------------------
        // Check for the optimum in last cell/column.
        // ----------------------------------------------------------------------------

        this->check_score_of_cells_in_last_column(alignment_column, this->alignment_state);
        this->check_score_of_last_cell(*alignment_column_it, this->alignment_state);

        // ----------------------------------------------------------------------------
        // Build the alignment result
        // ----------------------------------------------------------------------------

        static_assert(config_t::template exists<align_cfg::result>(),
                      "The configuration must contain an align_cfg::result element.");

        using result_t = typename align_result_selector<sequence1_t, sequence2_t, config_t>::type;
        result_t res{};

        res.id = idx;

        using result_config_t = std::remove_reference_t<
                decltype(seqan3::get<align_cfg::result>(std::declval<config_t>()).value)>;

        // Choose what needs to be computed.
        if constexpr (result_config_t::rank >= 0)  // compute score
            res.score = this->alignment_state.optimum.score;

        if constexpr (result_config_t::rank >= 1)  // compute back coordinate
        {
            res.back_coordinate = this->alignment_state.optimum.coordinate;
            // At some point this needs to be refactored so that it is not necessary to adapt the coordinate.
            if constexpr (is_banded)
                res.back_coordinate.second += res.back_coordinate.first - this->trace_matrix.band_col_index;
        }

        if constexpr (result_config_t::rank >= 2) // compute front coordinate
        {
            // Get a aligned sequence builder for banded or un-banded case.
            aligned_sequence_builder builder{sequence1, sequence2};

            auto trace_res = builder(this->trace_matrix.trace_path(this->alignment_state.optimum.coordinate));
            res.front_coordinate.first = trace_res.first_sequence_slice_positions.first;
            res.front_coordinate.second = trace_res.second_sequence_slice_positions.first;

            if constexpr (result_config_t::rank == 3) // compute alignment
                res.alignment = std::move(trace_res.alignment);
        }

        // Store the matrices in debug mode.
        if constexpr (is_debug_mode)
        {
            res.score_debug_matrix = std::move(score_debug_matrix);
            if constexpr (result_config_t::rank == 3) // compute alignment
                res.trace_debug_matrix = std::move(trace_debug_matrix);
        }

        return res;
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
        if constexpr (is_banded)
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
    //!\brief The scoring scheme used for this alignment algorithm.
    score_scheme_t score_scheme{};
    //!\brief Stores the currently processed alignment column.
    alignment_column_t alignment_column{};
    //!\brief Stores the state of the currently processed alignment column.
    alignment_column_iterator_t alignment_column_it{};
    //!\brief The debug matrix for the scores.
    score_debug_matrix_t score_debug_matrix{};
    //!\brief The debug matrix for the traces.
    trace_debug_matrix_t trace_debug_matrix{};
};

} // namespace seqan3::detail
