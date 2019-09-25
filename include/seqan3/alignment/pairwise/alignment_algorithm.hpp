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

#include <type_traits>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/exception.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alignment/matrix/detail/aligned_sequence_builder.hpp>
#include <seqan3/alignment/pairwise/detail/alignment_algorithm_cache.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_init_policy.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_policy.hpp>
#include <seqan3/alignment/pairwise/policy/unbanded_score_dp_matrix_policy.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>

#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/type_traits/deferred_crtp_base.hpp>
#include <seqan3/range/views/drop.hpp>
#include <seqan3/range/views/take_exactly.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/concepts>
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

    //!\brief The type of the score debug matrix type.
    using score_debug_matrix_t =
        std::conditional_t<is_debug_mode,
                           two_dimensional_matrix<std::optional<int32_t>,
                                                  std::allocator<std::optional<int32_t>>,
                                                  matrix_major_order::column>,
                           empty_type>;
    //!\brief The type of the trace debug matrix type.
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
     * Maintains a copy of the configuration object on the heap using a std::shared_ptr.
     */
    explicit constexpr alignment_algorithm(config_t const & cfg) : cfg_ptr{new config_t(cfg)}
    {}
    //!\}

    /*!\brief Invokes the actual alignment computation given two sequences.
     * \tparam    first_range_t  The type of the first sequence (or packed sequences); must model std::forward_range.
     * \tparam    second_range_t The type of the second sequence (or packed sequences); must model std::forward_range.
     *
     * \param[in] idx            The index of the current processed sequence pair.
     * \param[in] first_range    The first sequence (or packed sequences).
     * \param[in] second_range   The second sequence (or packed sequences).
     *
     * \returns A seqan3::alignment_result with the requested alignment outcomes.
     *
     * \details
     *
     * The algorithm always computes a pairwise alignment of two sequences over either a regular alphabet or
     * packed alphabets in a SIMD vector. In the latter case an inter-vectorisation layout
     * is used to compute l many pairwise alignments in parallel using special extended register instructions.
     *
     * ### Exception
     *
     * Strong exception guarantee.
     *
     * ### Thread-safety
     *
     * Calls to this functions in a concurrent environment are not thread safe. Instead use a copy of the alignment
     * algorithm type.
     *
     * ### Complexity
     *
     * The code always runs in \f$ O(N^2) \f$ time and depending on the configuration requires at least \f$ O(N) \f$
     * and at most \f$ O(N^2) \f$ space.
     */
    template <std::ranges::forward_range first_range_t, std::ranges::forward_range second_range_t>
    auto operator()(size_t const idx, first_range_t && first_range, second_range_t && second_range)
        requires !is_banded
    {
        assert(cfg_ptr != nullptr);

        if constexpr (is_debug_mode)
            initialise_debug_matrices(first_range, second_range);

        // ----------------------------------------------------------------------------
        // Initialise dp algorithm.
        // ----------------------------------------------------------------------------

        is_fst_range_empty = std::ranges::empty(first_range);
        is_snd_range_empty = std::ranges::empty(second_range);

        // We need to allocate the score_matrix and maybe the trace_matrix.
        this->allocate_matrix(first_range, second_range);

        // Get the selected or default gap scheme.
        auto cache = this->make_cache(
                cfg_ptr->template value_or<align_cfg::gap>(gap_scheme{gap_score{-1}, gap_open_score{-10}}));

        initialise_matrix(cache);

        // ----------------------------------------------------------------------------
        // Compute the unbanded alignment.
        // ----------------------------------------------------------------------------

        compute_matrix(first_range, second_range, cache);

        // ----------------------------------------------------------------------------
        // Make the alignment result.
        // ----------------------------------------------------------------------------

        return alignment_result{make_result(idx, first_range, second_range, cache)};
    }

    /*!\brief Invokes the banded alignment computation given two sequences.
     * \tparam    first_range_t  The type of the first sequence (or packed sequences); must model std::forward_range.
     * \tparam    second_range_t The type of the second sequence (or packed sequences); must model std::forward_range.
     *
     * \param[in] idx            The index of the current processed sequence pair.
     * \param[in] first_range    The first sequence (or packed sequences).
     * \param[in] second_range   The second sequence (or packed sequences).
     *
     * \returns A seqan3::alignment_result with the requested alignment outcomes.
     *
     * \details
     *
     * Computes the k-banded alignment. This function is only available if the alignment configuration was configured
     * with seqan3::align_cfg::band. The algorithm always computes a pairwise alignment of two sequences over either a
     * regular alphabet or packed alphabets in a SIMD vector. In the latter case an inter-vectorisation layout
     * is used to compute l many pairwise alignments in parallel using special extended register instructions.
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
     * The code always runs in \f$ O(N*k) \f$ time and depending on the configuration requires at least \f$ O(k) \f$
     * and at most \f$ O(N*k) \f$ space.
     */
    template <std::ranges::forward_range first_range_t, std::ranges::forward_range second_range_t>
    auto operator()(size_t const idx, first_range_t && first_range, second_range_t && second_range)
        requires is_banded
    {
        assert(cfg_ptr != nullptr);

        using std::get;

        static_assert(config_t::template exists<align_cfg::band>(),
                      "The band configuration is required for the banded alignment algorithm.");

        auto const & band = get<align_cfg::band>(*cfg_ptr).value;

        // ----------------------------------------------------------------------------
        // Check valid band settings.
        // ----------------------------------------------------------------------------

        using diff_type = typename std::iterator_traits<std::ranges::iterator_t<first_range_t>>::difference_type;
        static_assert(std::is_signed_v<diff_type>,  "Only signed types can be used to test the band parameters.");

        if (static_cast<diff_type>(band.lower_bound) > std::ranges::distance(first_range))
        {
            throw invalid_alignment_configuration
            {
                "Invalid band error: The lower bound excludes the whole alignment matrix."
            };
        }

        if (static_cast<diff_type>(band.upper_bound) < -std::ranges::distance(second_range))
        {
            throw invalid_alignment_configuration
            {
                "Invalid band error: The upper bound excludes the whole alignment matrix."
            };
        }

        if constexpr (is_debug_mode)
            initialise_debug_matrices(first_range, second_range);

        // ----------------------------------------------------------------------------
        // Initialise dp algorithm.
        // ----------------------------------------------------------------------------

        // Trim the sequences according to special band settings.
        auto [first_slice, second_slice] = this->slice_sequences(first_range, second_range, band);

        is_fst_range_empty = std::ranges::empty(first_slice);
        is_snd_range_empty = std::ranges::empty(second_slice);

        auto cache = this->make_cache(
            cfg_ptr->template value_or<align_cfg::gap>(gap_scheme{gap_score{-1}, gap_open_score{-10}}));
        // In the allocation we need to trim
        // But then the builder needs the sequences
        this->allocate_matrix(first_slice, second_slice, band, cache);

        initialise_matrix(cache);

        // ----------------------------------------------------------------------------
        // Compute the banded alignment.
        // ----------------------------------------------------------------------------

        compute_banded_matrix(first_slice, second_slice, cache);

        // ----------------------------------------------------------------------------
        // Cleanup and optionally compute the traceback.
        // ----------------------------------------------------------------------------

        return alignment_result{make_result(idx, first_slice, second_slice, cache)};
    }

private:

    //!\brief Initialises the debug matrices for the given sequences.
    template <typename first_range_t, typename second_range_t>
    constexpr void initialise_debug_matrices(first_range_t & first_range, second_range_t & second_range)
    {
        size_t rows = std::ranges::distance(second_range) + 1;
        size_t cols = std::ranges::distance(first_range) + 1;
        score_debug_matrix = score_debug_matrix_t{number_rows{rows}, number_cols{cols}};
        trace_debug_matrix = trace_debug_matrix_t{number_rows{rows}, number_cols{cols}};
    }

    /*!\brief Initialises the first column of the dynamic programming matrix.
     * \tparam         cache_t The cache type.
     * \param[in,out]  cache   The cache holding hot variables.
     */
    template <typename cache_t>
    void initialise_matrix(cache_t & cache)
    {
        // Get the initial column.
        auto zip_column = views::zip(*this->score_matrix_iter, *this->trace_matrix_iter);
        assert(!zip_column.empty()); // Must contain at least one element.

        // Initialise first cell.
        this->init_origin_cell(zip_column.front(), cache);
        // Take everything but last cell. Note the zip view does not preserve sized range at the moment.
        auto column = zip_column | views::take_exactly(std::ranges::size(*this->score_matrix_iter) - 1);
        if (column.empty()) // TODO [[unlikely]]
        {
            if (is_snd_range_empty)  // if second range is empty we are at the last row.
                this->check_score_last_row(zip_column.front(), cache);

            if (is_fst_range_empty) // If this was the only column to compute, track the last column/cell.
                this->check_score_last_column_or_cell(zip_column.front(), cache);
        }
        else
        {
            // Go to next cell and compute all inner cells.
            auto column_iter = std::ranges::next(column.begin());
            for (; column_iter != column.end(); ++column_iter)
                this->init_column_cell(*column_iter, cache);

            this->init_column_cell(*column_iter, cache);

            // Check if we need to track the last cell.
            if constexpr (is_banded)
            {
                if (static_cast<size_t>(this->score_matrix.band_row_index) == this->score_matrix.num_rows - 1) // TODO [[unlikely]]
                    this->check_score_last_row(*column_iter, cache);  // band row index is equal to the number of rows
            }
            else
            {
                this->check_score_last_row(*column_iter, cache);
            }

            if (is_fst_range_empty) // [[unlikely]] If this was the only column to compute, track the last column/cell.
                this->check_score_last_column_or_cell(*column_iter, cache);
        }

        if constexpr (is_debug_mode)
            dump_alignment_column(zip_column);
    }

    /*!\brief Compute the alignment by iterating over the dynamic programming matrix in a column wise manner.
     * \tparam        first_range_t  The type of the first sequence (or packed sequences).
     * \tparam        second_range_t The type of the second sequence (or packed sequences).
     * \tparam        cache_t        The type of the cache.
     * \param[in]     first_range    The first sequence.
     * \param[in]     second_range   The second sequence.
     * \param[in,out] cache          The cache holding hot variables.
     */
    template <typename first_range_t,
              typename second_range_t,
              typename cache_t>
    void compute_matrix(first_range_t & first_range,
                        second_range_t & second_range,
                        cache_t & cache)
    {
        using zip_column_t = decltype(views::zip(*this->score_matrix_iter, *this->trace_matrix_iter));
        using column_iter_t = std::ranges::iterator_t<decltype(std::declval<zip_column_t>() | views::take_exactly(1))>;

        if (is_fst_range_empty)  // Already tracked the max in initialisation.
            return;

        auto const & score_scheme = get<align_cfg::scoring>(*cfg_ptr).value;

        column_iter_t column_iter{};

        for (auto const seq1_value : first_range)
        {
            // Move to the next column in the alignment matrix and initialise first cell.
            zip_column_t zip_column = views::zip(*++this->score_matrix_iter, *++this->trace_matrix_iter);
            assert(!zip_column.empty());

            // If the column is empty (the second range is empty)
            if (is_snd_range_empty) // TODO [[unlikely]]
            {
                this->init_row_cell(zip_column.front(), cache);
                this->check_score_last_row(zip_column.front(), cache);
            }
            else
            {
                // Compute the column
                auto column = zip_column | views::take_exactly(std::ranges::size(*this->score_matrix_iter) - 1);
                assert(!column.empty()); // Must at least contain one element.
                column_iter = column.begin();

                // Initialise first row.
                this->init_row_cell(*column_iter, cache);

                ++column_iter;
                auto second_range_it = std::ranges::begin(second_range);
                // Compute all inner cells of the column.
                for (; column_iter != column.end(); ++column_iter, ++second_range_it)
                    this->compute_cell(*column_iter, cache, score_scheme.score(seq1_value, *second_range_it));

                assert(second_range_it != std::ranges::end(second_range));
                assert(second_range_it + 1 == std::ranges::end(second_range));
                // Compute last cell of column.
                this->compute_cell(*column_iter, cache, score_scheme.score(seq1_value, *second_range_it));
                this->check_score_last_row(*column_iter, cache);
            }

            if constexpr (is_debug_mode)
                dump_alignment_column(zip_column);
        }

        assert(this->score_matrix_iter != this->score_matrix.end());
        assert(this->trace_matrix_iter != this->trace_matrix.end());

        if (is_snd_range_empty) // TODO [[unlikely]]
        {
            this->check_score_last_column_or_cell(
                views::zip(*this->score_matrix_iter, *this->trace_matrix_iter).front(), cache);
        }
        else
        {
            this->check_score_last_column_or_cell(*column_iter, cache);
        }
    }

    /*!\brief Compute the alignment by iterating over the banded dynamic programming matrix in a column wise manner.
     * \tparam        first_range_t  The type of the first sequence (or packed sequences).
     * \tparam        second_range_t The type of the second sequence (or packed sequences).
     * \tparam        cache_t        The type of the cache.
     * \param[in]     first_range    The first sequence.
     * \param[in]     second_range   The second sequence.
     * \param[in,out] cache          The cache holding hot variables.
     */
    template <typename first_range_t,
              typename second_range_t,
              typename cache_t>
    void compute_banded_matrix(first_range_t & first_range,
                               second_range_t & second_range,
                               cache_t & cache)
    {
        using zip_column_t = decltype(views::zip(*this->score_matrix_iter, *this->trace_matrix_iter));
        using column_iter_t = std::ranges::iterator_t<decltype(std::declval<zip_column_t>() | views::take_exactly(1))>;

        if (is_fst_range_empty)  // Already tracked the max in initialisation.
            return;

        auto const & score_scheme = get<align_cfg::scoring>(*cfg_ptr).value;
        auto last_row = std::ranges::next(std::ranges::begin(second_range), std::ranges::size(second_range) - 1);

        // ----------------------------------------------------------------------------
        // 1st phase: Iterate as long as the band intersects with the first row.
        // ----------------------------------------------------------------------------

        column_iter_t column_iter{};
        for (auto const first_range_value : first_range | views::take_exactly(this->score_matrix.band_col_index))
        {
            // Prepare next column
            zip_column_t zip_column = views::zip(*++this->score_matrix_iter, *++this->trace_matrix_iter);
            assert(!zip_column.empty());

            auto column = zip_column | views::take_exactly(std::ranges::size(*(this->score_matrix_iter)) - 1);
            if (column.empty()) // TODO [[unlikely]]
            {
                this->init_row_cell(zip_column.front(), cache);
                this->check_score_last_row(zip_column.front(), cache);
            }
            else
            {
                column_iter = column.begin();
                // Initialise first row.
                this->init_row_cell(*column_iter, cache);
                ++column_iter;
                // Compute all inner cells of the column.
                auto second_range_it = std::ranges::begin(second_range);
                for (; column_iter != column.end(); ++column_iter, ++second_range_it)
                    this->compute_cell(*column_iter, cache, score_scheme.score(first_range_value, *second_range_it));

                // Compute last cell of column.
                assert(second_range_it != std::ranges::end(second_range));
                this->compute_cell(*column_iter, cache, score_scheme.score(first_range_value, *second_range_it));

                // Check if we reached the last row.
                if (second_range_it == last_row)  // TODO [[unlikely]]
                    this->check_score_last_row(*column_iter, cache);
            }

            if constexpr (is_debug_mode)
                dump_alignment_column(zip_column);
        }

        // ----------------------------------------------------------------------------
        // 2nd phase: Iterate until the end of the matrix.
        // ----------------------------------------------------------------------------

        bool last_column_has_one_cell = false;  // In case last column has only one cell.
        auto begin_row = std::ranges::begin(second_range);
        for (auto const first_range_value : first_range | views::drop(this->score_matrix.band_col_index))
        {
            // Prepare next column.
            zip_column_t zip_column = views::zip(*++this->score_matrix_iter, *++this->trace_matrix_iter);
            assert(!zip_column.empty());

            auto column = zip_column | views::take_exactly(std::ranges::size(*(this->score_matrix_iter)) - 1);
            if (column.empty()) // TODO [[unlikely]]
            { // Either last column at end of band has size one or band has size 1.
                this->compute_first_band_cell(zip_column.front(),
                                              cache,
                                              score_scheme.score(first_range_value, *begin_row));
                last_column_has_one_cell = true;
            }
            else
            {
                column_iter = column.begin();

                // Compute first cell of band inside the matrix
                this->compute_first_band_cell(*column_iter,
                                              cache,
                                              score_scheme.score(first_range_value, *begin_row));
                ++column_iter;

                // Compute all inner cells of the column.
                auto second_range_it = std::ranges::next(begin_row);
                for (; column_iter != column.end(); ++column_iter, ++second_range_it)
                    this->compute_cell(*column_iter, cache, score_scheme.score(first_range_value, *second_range_it));

                assert(second_range_it != std::ranges::end(second_range));
                // Compute last cell of column.
                this->compute_cell(*column_iter, cache, score_scheme.score(first_range_value, *second_range_it));

                // Check if last row of alignment matrix was reached.
                if (second_range_it == last_row)  // TODO [[unlikely]]
                    this->check_score_last_row(*column_iter, cache);
            }
            ++begin_row;  // Move down the band and start at next row.

            if constexpr (is_debug_mode)
                dump_alignment_column(zip_column);
        }

        if (last_column_has_one_cell) // TODO [[unlikely]]
        {
            this->check_score_last_column_or_cell(
                    views::zip(*this->score_matrix_iter, *this->trace_matrix_iter).front(), cache);
        }
        else
        {
            this->check_score_last_column_or_cell(*column_iter, cache);
        }
    }

    /*!\brief Creates a new alignment result from the current alignment optimum and for the given pair of ranges.
     * \param[in] index The internal index used for this pair of sequences.
     * \param[in] fst_range The first range to get the alignment for if requested.
     * \param[in] snd_range The second range to get the alignment for if requested.
     * \param[in] cache The current alignment cache containing the computed optimum.
     *
     * \returns A seqan3::alignment_result with the requested alignment outcomes.
     *
     * \details
     *
     * Stores the computed score and id in a seqan3::alignment_result object and optionally stores the
     * subrange positions of the alignment and/or computes the alignment based on the trace matrix.
     */
    template <typename index_t, typename fst_range_t, typename snd_range_t, typename cache_t>
    auto make_result(index_t index, fst_range_t & fst_range, snd_range_t & snd_range, cache_t const & cache)
    {
        static_assert(config_t::template exists<align_cfg::result>(),
                "The configuration must contain an align_cfg::result element.");

        using result_t = typename align_result_selector<fst_range_t, snd_range_t, config_t>::type;
        result_t res{};

        res.id = index;

        using result_config_t = std::remove_reference_t<
                decltype(seqan3::get<align_cfg::result>(std::declval<config_t>()).value)>;

        // Choose what needs to be computed.
        if constexpr (result_config_t::rank >= 0)  // compute score
            res.score = cache.optimum.score;

        if constexpr (result_config_t::rank >= 1)  // compute back coordinate
        {
            res.back_coordinate = cache.optimum.coordinate;
            // At some point this needs to be refactored so that it is not necessary to adapt the coordinate.
            if constexpr (is_banded)
                res.back_coordinate.second += res.back_coordinate.first - this->trace_matrix.band_col_index;
        }

        if constexpr (result_config_t::rank >= 2) // compute front coordinate
        {
            // Get a aligned sequence builder for banded or un-banded case.
            aligned_sequence_builder builder{fst_range, snd_range};

            auto trace_res = builder(this->trace_matrix.trace_path(cache.optimum.coordinate));
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
     * \param[in] column The current column.
     *
     * \details
     *
     * Copies the current score and if configured the current trace column into the local debug score and trace matrix.
     * In the banded case the full matrix will be allocated with std::optional values and only the respective parts
     * of the matrix that are present in the band are filled.
     */
    template <typename zipped_column_t>
    void dump_alignment_column(zipped_column_t column)
    {
        using std::get;
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
    //!\brief The debug matrix for the scores.
    score_debug_matrix_t score_debug_matrix{};
    //!\brief The debug matrix for the traces.
    trace_debug_matrix_t trace_debug_matrix{};
    //!\brief Flag to check if first range is empty.
    bool is_fst_range_empty{false};
    //!\brief Flag to check if second range is empty.
    bool is_snd_range_empty{false};
};

} // namespace seqan3::detail
