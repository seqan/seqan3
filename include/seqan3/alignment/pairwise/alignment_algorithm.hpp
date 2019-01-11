// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_algorithm.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/drop_exactly.hpp>
#include <range/v3/view/zip.hpp>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/exception.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_init_policy.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_policy.hpp>
#include <seqan3/alignment/pairwise/policy/unbanded_dp_matrix_policy.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>

#include <seqan3/core/metafunction/deferred_crtp_base.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/range/view/get.hpp>
#include <seqan3/range/view/take_exactly.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief The alignment algorithm type to compute standard pairwise alignment using dynamic programming.
 * \implements std::Invocable
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
    static constexpr bool is_banded = std::remove_reference_t<config_t>::template exists<align_cfg::band>();

public:
    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr alignment_algorithm()                                        = default;
    constexpr alignment_algorithm(alignment_algorithm const &)             = default;
    constexpr alignment_algorithm(alignment_algorithm &&)                  = default;
    constexpr alignment_algorithm & operator=(alignment_algorithm const &) = default;
    constexpr alignment_algorithm & operator=(alignment_algorithm &&)      = default;
    ~alignment_algorithm()                                                 = default;

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
     * \tparam    first_range_t  The type of the first sequence (or packed sequences); must model std::ForwardRange.
     * \tparam    second_range_t The type of the second sequence (or packed sequences); must model std::ForwardRange.
     * \param[in] first_range    The first sequence (or packed sequences).
     * \param[in] second_range   The second sequence (or packed sequences).
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
    template <std::ranges::ForwardRange first_range_t, std::ranges::ForwardRange second_range_t>
    auto operator()(first_range_t const & first_range, second_range_t const & second_range)
        requires !is_banded
    {
        assert(cfg_ptr != nullptr);

        // We need to allocate the score_matrix and maybe the trace_matrix.
        this->allocate_matrix(first_range, second_range);

        // Initialise cache variables to keep frequently used variables close to the CPU registers.
        auto cache = this->make_cache(cfg_ptr->template value_or<align_cfg::gap>(gap_scheme{gap_score{-1},
                                                                                            gap_open_score{-10}}));

        initialise_matrix(cache);

        compute_matrix(first_range, second_range, cache);

        using result_t = typename align_result_selector<first_range_t, second_range_t, config_t>::type;
        result_t res{};

        // Choose what needs to be computed.
        if constexpr (config_t::template exists<align_cfg::result<detail::with_score_type>>())
        {
            res.score = get<3>(cache);
        }
        return align_result{res};
    }

    /*!\brief Invokes the banded alignment computation given two sequences.
     * \tparam    first_range_t  The type of the first sequence (or packed sequences); must model std::ForwardRange.
     * \tparam    second_range_t The type of the second sequence (or packed sequences); must model std::ForwardRange.
     * \param[in] first_range    The first sequence (or packed sequences).
     * \param[in] second_range   The second sequence (or packed sequences).
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
    template <std::ranges::ForwardRange first_range_t, std::ranges::ForwardRange second_range_t>
    auto operator()(first_range_t const & first_range, second_range_t const & second_range)
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

        if (static_cast<diff_type>(band.lower_bound) >
            std::ranges::distance(seqan3::begin(first_range), seqan3::end(first_range)))
        {
            throw invalid_alignment_configuration
            {
                "The lower bound is bigger than the size of the first sequence."
            };
        }

        if (static_cast<diff_type>(std::abs(band.upper_bound)) >
            std::ranges::distance(seqan3::begin(second_range), seqan3::end(second_range)))
        {
            throw invalid_alignment_configuration
            {
                "The upper bound is smaller than the size of the second sequence."
            };
        }

        // ----------------------------------------------------------------------------
        // Initialise dp algorithm.
        // ----------------------------------------------------------------------------

        // Trim the sequences according to special band settings.
        auto [trimmed_first_range, trimmed_second_range] =
            this->trim_sequences(first_range, second_range, band);

        this->allocate_matrix(trimmed_first_range, trimmed_second_range, band);

        // Use default gap if not set from outside.
        auto const & gap = cfg_ptr->template value_or<align_cfg::gap>(gap_scheme{gap_score{-1}, gap_open_score{-10}});

        // Initialise cache variables to keep frequently used variables close to the CPU registers.
        auto cache = this->make_cache(gap);

        initialise_matrix(cache);

        // ----------------------------------------------------------------------------
        // Compute the banded alignment.
        // ----------------------------------------------------------------------------

        compute_banded_matrix(trimmed_first_range, trimmed_second_range, cache);

        // ----------------------------------------------------------------------------
        // Cleanup and optionally compute the traceback.
        // ----------------------------------------------------------------------------

        using result_t = typename align_result_selector<first_range_t, second_range_t, config_t>::type;
        result_t res{};

        // Balance the score with possible leading/trailing gaps depending on the
        // band settings.
        this->balance_leading_gaps(get<3>(cache), band, gap);

        this->balance_trailing_gaps(get<3>(cache),
                                    this->dimension_first_range,
                                    this->dimension_second_range,
                                    band,
                                    gap);

        if constexpr (config_t::template exists<align_cfg::result<detail::with_score_type>>())
        {
            res.score = get<3>(cache);
        }
        return align_result{res};
    }
private:

    /*!\brief Initialises the first column of the dynamic programming matrix.
     * \tparam         cache_t The cache type.
     * \param[in,out]  cache   The cache holding hot variables.
     */
    template <typename cache_t>
    void initialise_matrix(cache_t & cache)
    {
        auto col = this->current_column();
        auto score_col = col | view_get_score_column;

        this->init_origin_cell(*seqan3::begin(col), cache);

        ranges::for_each(col | ranges::view::drop_exactly(1), [&cache, this](auto && cell)
        {
            this->init_column_cell(std::forward<decltype(cell)>(cell), cache);
        });

        if constexpr (is_banded)
            this->check_score_last_row(get<0>(get<0>(*std::ranges::prev(seqan3::end(score_col)))), get<3>(cache));
        else
            this->check_score_last_row(get<0>(*std::ranges::prev(seqan3::end(score_col))), get<3>(cache));
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
    void compute_matrix(first_range_t const & first_range,
                        second_range_t const & second_range,
                        cache_t & cache)
    {
        using std::get;
        auto const & score_scheme = get<align_cfg::scoring>(*cfg_ptr).value;
        ranges::for_each(first_range, [&, this](auto seq1_value)
        {
            auto col = this->current_column();
            auto score_col = col | view_get_score_column;
            this->init_row_cell(*seqan3::begin(col), cache);

            auto second_range_it = std::ranges::begin(second_range);
            ranges::for_each(col | ranges::view::drop_exactly(1), [&, this] (auto && cell)
            {
                this->compute_cell(cell, cache, score_scheme.score(seq1_value, *second_range_it));
                ++second_range_it;
            });
            this->check_score_last_row(get<0>(*std::ranges::prev(seqan3::end(score_col))), get<3>(cache));
        });
        this->check_score_last_column(this->current_column() | view_get_score_column, get<3>(cache));
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
    void compute_banded_matrix(first_range_t const & first_range,
                               second_range_t const & second_range,
                               cache_t & cache)
    {
        auto const & score_scheme = get<align_cfg::scoring>(*cfg_ptr).value;
        // ----------------------------------------------------------------------------
        // 1st phase: Iterate as long as the band intersects with the first row.
        // ----------------------------------------------------------------------------

        ranges::for_each(first_range | view::take_exactly(this->band_column_index), [&, this](auto first_range_value)
        {
            this->next_column(); // Move to the next column.
            auto col = this->current_column();
            auto score_col = col | view_get_score_column;

            this->init_row_cell(*seqan3::begin(col), cache); // initialise first row of dp matrix.

            auto second_range_it = seqan3::begin(second_range);
            ranges::for_each(col | ranges::view::drop_exactly(1), [&, this] (auto && cell)
            {
                this->compute_cell(std::forward<decltype(cell)>(cell),
                                   cache,
                                   score_scheme.score(first_range_value, *second_range_it));
                ++second_range_it;
            });

            if (this->band_touches_last_row())  // TODO [[unlikely]]
                this->check_score_last_row(get<0>(get<0>(*std::ranges::prev(seqan3::end(score_col)))), get<3>(cache));
        });

        // ----------------------------------------------------------------------------
        // 2nd phase: Iterate until the end of the matrix.
        // ----------------------------------------------------------------------------

        // Drop the first columns from the 1st phase.
        ranges::for_each(first_range | ranges::view::drop_exactly(this->band_column_index),
        [&, this](auto first_range_value)
        {
            this->next_column(); // Move to the next column.
            auto col = this->current_column();
            auto score_col = col | view_get_score_column;

            // Move the second_range_it to the correct position depending on the current band position.
            auto second_range_it = std::ranges::begin(second_range);
            std::advance(second_range_it, this->second_range_begin_offset());

            // Initialise the first cell of the current band.
            this->compute_first_band_cell(*std::ranges::begin(col),
                                          cache,
                                          score_scheme.score(first_range_value, *second_range_it));

            // Process rest of current column band.
            ++second_range_it;
            ranges::for_each(col | ranges::view::drop_exactly(1), [&, this] (auto && cell)
            {
                this->compute_cell(std::forward<decltype(cell)>(cell),
                                   cache,
                                   score_scheme.score(first_range_value, *second_range_it));
                ++second_range_it;
            });

            if (this->band_touches_last_row()) // TODO [[unlikely]]
                this->check_score_last_row(get<0>(get<0>(*std::ranges::prev(seqan3::end(score_col)))), get<3>(cache));
        });
        // We need to call view_get_score_column twice. The first time we access the score column which is a
        // zipped range in the banded case and the second time to call the actual score column.
        this->check_score_last_column(this->current_column() | view_get_score_column | view_get_score_column,
                                      get<3>(cache));
    }

    //!\brief The alignment configuration stored on the heap.
    std::shared_ptr<config_t> cfg_ptr{};
};

} // namespace seqan3::detail
