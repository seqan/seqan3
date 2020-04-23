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
 *
 * \details
 *
 * ### Configuration
 *
 * The first template argument is the type of the alignment configuration which was used to configure the alignment
 * algorithm type within the seqan3::detail::alignment_configurator. The algorithm computes a column based
 * dynamic programming matrix given two sequences. After the computation a user defined callback function is invoked
 * with the computed seqan3::alignment_result.
 */
template <typename alignment_configuration_t>
//!\cond
    requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
//!\endcond
class pairwise_alignment_algorithm
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
    //!\brief The configured gap scheme.
    gap_scheme<int8_t> m_gap_scheme{};

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
    pairwise_alignment_algorithm(alignment_configuration_t const & config)
    {
        m_scoring_scheme = seqan3::get<align_cfg::scoring>(config).value;
        m_gap_scheme = config.template value_or<align_cfg::gap>(gap_scheme<int8_t>{gap_score{-1}, gap_open_score{-10}});
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
    /*!\brief Compute the actual alignment by given the two sequences.
     * \tparam sequence1_t The type of the first sequence; must model std::ranges::forward_range.
     * \tparam sequence2_t The type of the second sequence; std::ranges::forward_range.
     *
     * \param[in] sequence1 The first sequence.
     * \param[in] sequence2 The second sequence.
     */
    template <std::ranges::forward_range sequence1_t, std::ranges::forward_range sequence2_t>
    int32_t compute_matrix(sequence1_t & sequence1, sequence2_t & sequence2)
    {
        // Use gap costs locally.
        int32_t gap_extension{m_gap_scheme.get_gap_score()};
        int32_t gap_open{m_gap_scheme.get_gap_open_score() + gap_extension};

        std::vector<int32_t> optimal_column{};
        std::vector<int32_t> horizontal_column{};

        // Initialise matrix.
        optimal_column.clear();
        horizontal_column.clear();
        optimal_column.resize(sequence2.size() + 1, 0);
        horizontal_column.resize(sequence2.size() + 1, 0);
        int32_t diagonal{};
        int32_t vertical{gap_open};

        // Initialise the first column.
        for (auto && [opt, hor] : views::zip(optimal_column, horizontal_column) | views::drop(1))
        {
            opt = vertical;
            hor = opt + gap_open;
            vertical += gap_extension;
        }

        // Compute the matrix
        for (auto it_col = sequence1.begin(); it_col != sequence1.end(); ++it_col)
        {
            // Initialise first cell of optimal_column.
            auto opt_it = optimal_column.begin();
            auto hor_it = horizontal_column.begin();

            diagonal = *opt_it;  // Cache the diagonal for next cell.
            *opt_it =  gap_open + gap_extension * (it_col - sequence1.begin()); // Initialise the main score.
            *hor_it = *opt_it; // Initialise the horizontal score.
            vertical = *opt_it + gap_open; // Initialise the vertical score.
            // Go to next cell.
            ++opt_it;
            ++hor_it;
            for (auto it_row = sequence2.begin(); it_row != sequence2.end(); ++it_row, ++opt_it, ++hor_it)
            {
                // Precompute the diagonal score.
                int32_t tmp = diagonal + m_scoring_scheme.score(*it_col, *it_row);

                tmp = (tmp < vertical) ? vertical : tmp;
                tmp = (tmp < *hor_it) ? *hor_it : tmp;

                // Store the current max score.
                diagonal = *opt_it; // Cache the next diagonal before writing it.
                *opt_it = tmp; // Store the temporary result.

                tmp += gap_open;  // Add gap open costs.
                vertical += gap_extension;
                *hor_it += gap_extension;

                // Store the vertical and horizontal value in the next path.
                vertical = (vertical < tmp) ? tmp : vertical;
                *hor_it = (*hor_it < tmp) ? tmp : *hor_it;
            }
        }
        return optimal_column.back();
    }
};
} // namespace seqan3::detail
