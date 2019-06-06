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

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/drop_exactly.hpp>
#include <range/v3/view/zip.hpp>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/exception.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_init_policy.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_policy.hpp>
#include <seqan3/alignment/pairwise/policy/unbanded_score_dp_matrix_policy.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>

#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/type_traits/deferred_crtp_base.hpp>
#include <seqan3/range/view/get.hpp>
#include <seqan3/range/view/slice.hpp>
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
     * \tparam    first_range_t  The type of the first sequence (or packed sequences); must model std::ForwardRange.
     * \tparam    second_range_t The type of the second sequence (or packed sequences); must model std::ForwardRange.
     * \param[in] idx            The index of the current processed sequence pair.
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
    auto operator()(size_t const idx, first_range_t && first_range, second_range_t && second_range)
        requires !is_banded
    {
        assert(cfg_ptr != nullptr);

        // ----------------------------------------------------------------------------
        // Initialise dp algorithm.
        // ----------------------------------------------------------------------------

        // We need to allocate the score_matrix and maybe the trace_matrix.
        this->allocate_matrix(first_range, second_range);

        // Initialise cache variables to keep frequently used variables close to the CPU registers.
        auto cache = this->make_cache(cfg_ptr->template value_or<align_cfg::gap>(gap_scheme{gap_score{-1},
                                                                                            gap_open_score{-10}}));

        initialise_matrix(cache);

        // ----------------------------------------------------------------------------
        // Compute the unbanded alignment.
        // ----------------------------------------------------------------------------

        compute_matrix(first_range, second_range, cache);

        // ----------------------------------------------------------------------------
        // Cleanup and prepare the alignment result.
        // ----------------------------------------------------------------------------

        using result_t = typename align_result_selector<first_range_t, second_range_t, config_t>::type;
        result_t res{};

        res.id = idx;
        // Choose what needs to be computed.
        if constexpr (config_t::template exists<align_cfg::result<with_score_type>>())
        {
            res.score = get<3>(cache).score;
        }
        if constexpr (config_t::template exists<align_cfg::result<with_back_coordinate_type>>())
        {
            res.score = get<3>(cache).score;
            res.back_coordinate = get<3>(cache).coordinate;
        }
        if constexpr (config_t::template exists<align_cfg::result<with_front_coordinate_type>>())
        { // At the moment we also compute the traceback even if only the front coordinate was requested.
          // This can be later optimised by computing the reverse alignment in a narrow band in linear memory.
          // Especially for the SIMD version this might be more efficient.
            res.score = get<3>(cache).score;
            res.back_coordinate = get<3>(cache).coordinate;
            res.front_coordinate = get<0>(compute_traceback(first_range,
                                                            second_range,
                                                            get<3>(cache).coordinate));
        }
        if constexpr (config_t::template exists<align_cfg::result<with_alignment_type>>())
        {
            res.score = get<3>(cache).score;
            res.back_coordinate = get<3>(cache).coordinate;
            std::tie(res.front_coordinate, res.alignment) = compute_traceback(first_range,
                                                                              second_range,
                                                                              get<3>(cache).coordinate);
        }
        return alignment_result{res};
    }

    /*!\brief Invokes the banded alignment computation given two sequences.
     * \tparam    first_range_t  The type of the first sequence (or packed sequences); must model std::ForwardRange.
     * \tparam    second_range_t The type of the second sequence (or packed sequences); must model std::ForwardRange.
     * \param[in] idx            The index of the current processed sequence pair.
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

        res.id = idx;
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
            res.score = get<3>(cache).score;
        }
        if constexpr (config_t::template exists<align_cfg::result<with_back_coordinate_type>>())
        {
            res.score = get<3>(cache).score;
            res.back_coordinate = this->map_banded_coordinate_to_range_position(get<3>(cache).coordinate);
        }
        if constexpr (config_t::template exists<align_cfg::result<with_front_coordinate_type>>())
        { // At the moment we also compute the traceback even if only the front coordinate was requested.
          // This can be later optimised by computing the reverse alignment in linear memory from the maximum.
          // Especially for the SIMD version this might be more efficient.
            res.score = get<3>(cache).score;
            res.back_coordinate = this->map_banded_coordinate_to_range_position(get<3>(cache).coordinate);
            res.front_coordinate =
                get<0>(compute_traceback(first_range,
                                         second_range,
                                         get<3>(cache).coordinate));
        }
        if constexpr (config_t::template exists<align_cfg::result<with_alignment_type>>())
        {
            res.score = get<3>(cache).score;
            res.back_coordinate = this->map_banded_coordinate_to_range_position(get<3>(cache).coordinate);
            std::tie(res.front_coordinate, res.alignment) = compute_traceback(first_range,
                                                                              second_range,
                                                                              get<3>(cache).coordinate);
        }
        return alignment_result{res};
    }
private:

    /*!\brief Initialises the first column of the dynamic programming matrix.
     * \tparam         cache_t The cache type.
     * \param[in,out]  cache   The cache holding hot variables.
     */
    template <typename cache_t>
    void initialise_matrix(cache_t & cache)
    {
        // Get the current dynamic programming matrix.
        auto col = this->current_column();

        this->init_origin_cell(*std::ranges::begin(col), cache);

        ranges::for_each(col | ranges::view::drop_exactly(1), [&cache, this](auto && cell)
        {
            this->init_column_cell(std::forward<decltype(cell)>(cell), cache);
        });

        auto [cell, coordinate, trace] = *std::ranges::prev(std::ranges::end(col));
        (void) trace;
        if constexpr (is_banded)
        {
            alignment_optimum current{get<0>(get<0>(std::move(cell))), static_cast<alignment_coordinate>(coordinate)};
            this->check_score_last_row(current, get<3>(cache));
        }
        else
        {
            alignment_optimum current{get<0>(std::move(cell)), static_cast<alignment_coordinate>(coordinate)};
            this->check_score_last_row(current, get<3>(cache));
        }
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
        using std::get;
        auto const & score_scheme = get<align_cfg::scoring>(*cfg_ptr).value;
        ranges::for_each(first_range, [&, this](auto seq1_value)
        {
            // Move internal matrix to next column.
            this->go_next_column();

            auto col = this->current_column();
            this->init_row_cell(*std::ranges::begin(col), cache);

            auto second_range_it = std::ranges::begin(second_range);
            ranges::for_each(col | ranges::view::drop_exactly(1), [&, this] (auto && cell)
            {
                this->compute_cell(cell, cache, score_scheme.score(seq1_value, *second_range_it));
                ++second_range_it;
            });

            // Prepare last cell for tracking the optimum.
            auto [cell, coordinate, trace] = *std::ranges::prev(std::ranges::end(col));
            (void) trace;
            alignment_optimum current{get<0>(std::move(cell)), static_cast<alignment_coordinate>(coordinate)};
            this->check_score_last_row(current, get<3>(cache));
        });

        // Prepare the last column for tracking the optimum: Only get the current score cell and the coordinate.
        auto last_column_view = this->current_column() | std::view::transform([](auto && entry)
            {
            using std::get;
            return std::tuple{get<0>(std::forward<decltype(entry)>(entry)),
                              get<1>(std::forward<decltype(entry)>(entry))};
        });
        this->check_score_last_column(last_column_view, get<3>(cache));
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
        auto const & score_scheme = get<align_cfg::scoring>(*cfg_ptr).value;
        // ----------------------------------------------------------------------------
        // 1st phase: Iterate as long as the band intersects with the first row.
        // ----------------------------------------------------------------------------

        ranges::for_each(first_range | view::take_exactly(this->band_column_index), [&, this](auto first_range_value)
        {
            this->go_next_column(); // Move to the next column.
            auto col = this->current_column();
            this->init_row_cell(*std::ranges::begin(col), cache); // initialise first row of dp matrix.

            auto second_range_it = std::ranges::begin(second_range);
            ranges::for_each(col | ranges::view::drop_exactly(1), [&, this] (auto && cell)
            {
                this->compute_cell(std::forward<decltype(cell)>(cell),
                                   cache,
                                   score_scheme.score(first_range_value, *second_range_it));
                ++second_range_it;
            });

            if (this->band_touches_last_row())  // TODO [[unlikely]]
            {
                auto [cell, coordinate, trace] = *std::ranges::prev(std::ranges::end(col));
                (void) trace;
                alignment_optimum current{get<0>(get<0>(std::move(cell))),
                                          static_cast<alignment_coordinate>(coordinate)};
                this->check_score_last_row(current, get<3>(cache));
            }
        });

        // ----------------------------------------------------------------------------
        // 2nd phase: Iterate until the end of the matrix.
        // ----------------------------------------------------------------------------

        // Drop the first columns from the 1st phase.
        ranges::for_each(first_range | ranges::view::drop_exactly(this->band_column_index),
        [&, this](auto first_range_value)
        {
            this->go_next_column(); // Move to the next column.
            auto col = this->current_column();

            // Move the second_range_it to the correct position depending on the current band position.
            auto second_range_it = std::ranges::begin(second_range);
            std::ranges::advance(second_range_it, this->second_range_begin_offset());

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
            {
                auto [cell, coordinate, trace] = *std::ranges::prev(std::ranges::end(col));
                (void) trace;
                alignment_optimum current{get<0>(get<0>(std::move(cell))),
                                          static_cast<alignment_coordinate>(coordinate)};
                this->check_score_last_row(current, get<3>(cache));
            }
        });
        // Prepare the last column for tracking the optimum: Only get the current score cell and the coordinate.
        auto last_column_view = this->current_column() | std::view::transform([](auto && entry) {
            using std::get;
            return std::tuple{get<0>(get<0>(std::forward<decltype(entry)>(entry))),
                              get<1>(std::forward<decltype(entry)>(entry))};
        });
        this->check_score_last_column(last_column_view, get<3>(cache));
    }

    /*!\brief Computes the traceback if requested.
    * \tparam    first_range_t  The type of the first sequence (or packed sequences).
    * \tparam    second_range_t The type of the second sequence (or packed sequences).
    * \param[in] first_range    The first sequence.
    * \param[in] second_range   The second sequence.
    * \param[in] back_coordinate The back coordinate within the matrix where the traceback starts.
    *
    * \details
    *
    * First parses the traceback and computes the gap segments for the sequences. Then applies the gap segments
    * to the infix of the corresponding range and return the aligned sequence.
    */
    template <typename first_range_t, typename second_range_t>
    auto compute_traceback(first_range_t & first_range,
                           second_range_t & second_range,
                           alignment_coordinate back_coordinate)
    {
        // Parse the traceback
        auto [front_coordinate, first_gap_segments, second_gap_segments] = this->parse_traceback(back_coordinate);

        using result_type = typename align_result_selector<first_range_t, second_range_t, config_t>::type;
        using aligned_seq_type = decltype(result_type{}.alignment);

        // If we compute the traceback the aligned sequences must be provided.
        if constexpr (seqan3::TupleLike<aligned_seq_type>)
        {
            auto fill_aligned_sequence = [] (auto & aligned_sequence, auto & gap_segments, size_t const normalise)
            {
                assert(std::ranges::empty(gap_segments) || normalise <= gap_segments[0].position);

                size_t offset = 0;
                for (auto const & gap_elem : gap_segments)
                {
                    auto it = std::ranges::begin(aligned_sequence);
                    std::ranges::advance(it, (gap_elem.position - normalise) + offset);
                    insert_gap(aligned_sequence, it, gap_elem.size);
                    offset += gap_elem.size;
                }
            };

            // In banded case we need to refine the back coordinate to map to the correct position within the
            // second range.
            if constexpr (is_banded)
                back_coordinate = this->map_banded_coordinate_to_range_position(back_coordinate);

            // Get the subrange over the first sequence according to the front and back coordinate.
            auto first_subrange = view::slice(first_range, front_coordinate.first, back_coordinate.first);

            // Get the subrange over the second sequence according to the front and back coordinate.
            auto second_subrange = view::slice(second_range, front_coordinate.second, back_coordinate.second);

            // Create and fill the aligned_sequences.
            aligned_seq_type aligned_seq;
            assign_unaligned(std::get<0>(aligned_seq), first_subrange);
            assign_unaligned(std::get<1>(aligned_seq), second_subrange);
            fill_aligned_sequence(std::get<0>(aligned_seq), first_gap_segments, front_coordinate.first);
            fill_aligned_sequence(std::get<1>(aligned_seq), second_gap_segments, front_coordinate.second);

            return std::tuple{front_coordinate, aligned_seq};
        }
        else
        {
            return std::tuple{front_coordinate, std::ignore};
        }
    }

    //!\brief The alignment configuration stored on the heap.
    std::shared_ptr<config_t> cfg_ptr{};
};

} // namespace seqan3::detail
