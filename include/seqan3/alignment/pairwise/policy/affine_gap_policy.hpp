// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::affine_gap_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <limits>
#include <tuple>

#include <seqan3/alignment/matrix/alignment_optimum.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alignment/pairwise/detail/alignment_algorithm_cache.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// affine_gap_policy
// ----------------------------------------------------------------------------

/*!\brief The hot kernel implementation using affine gaps.
 * \ingroup alignment_policy
 * \tparam derived_t The derived alignment algorithm.
 * \tparam score_t The score type of the dynamic programming matrix.
 */
template <typename derived_t, typename score_t, typename align_local_t = std::false_type>
class affine_gap_policy
{
private:

    //!\brief Befriends the derived class to grant it access to the private members.
    friend derived_t;

    /*!\name Constructors, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr affine_gap_policy() noexcept = default;                                      //!< Defaulted
    constexpr affine_gap_policy(affine_gap_policy const &) noexcept = default;             //!< Defaulted
    constexpr affine_gap_policy(affine_gap_policy &&) noexcept = default;                  //!< Defaulted
    constexpr affine_gap_policy & operator=(affine_gap_policy const &) noexcept = default; //!< Defaulted
    constexpr affine_gap_policy & operator=(affine_gap_policy &&) noexcept = default;      //!< Defaulted
    ~affine_gap_policy() noexcept = default;                                               //!< Defaulted
    //!\}

    /*!\brief Computes the score of the current cell.
     * \tparam  cell_t     The type of the current cell.
     * \tparam  cache_type The type of the cache.
     * \param[in,out] current_cell The current cell in the dynamic programming matrix.
     * \param[in,out] cache        The cache storing hot helper variables.
     * \param[in]     score        The score of comparing the respective letters of the first and the second sequence.
     */
    template <typename cell_t, typename cache_t>
    constexpr void compute_cell(cell_t && current_cell,
                                cache_t & cache,
                                score_t const score) const noexcept
    {
        auto [score_cell, trace_cell] = current_cell;
        constexpr bool with_trace = !decays_to_ignore_v<std::remove_reference_t<decltype(trace_cell.current)>>;
        // Precompute the diagonal score.
        score_t tmp = score_cell.diagonal + score;
        // TODO: Differentiate between trace and no trace.
        if constexpr (with_trace)
        {
            tmp = (tmp < score_cell.up) ? (trace_cell.current = trace_cell.up, score_cell.up)
                                        : (trace_cell.current = trace_directions::diagonal | trace_cell.up, tmp);
            tmp = (tmp < score_cell.r_left) ? (trace_cell.current = trace_cell.r_left, score_cell.r_left)
                                            : (trace_cell.current |= trace_cell.r_left, tmp);
        }
        else
        {
            tmp = (tmp < score_cell.up) ? score_cell.up : tmp;
            tmp = (tmp < score_cell.r_left) ? score_cell.r_left : tmp;
        }

        if constexpr (align_local_t::value)
            tmp = (tmp < 0) ? (trace_cell.current = trace_directions::none, 0) : tmp;

        // Store the current max score.
        score_cell.current = tmp;
        // Check if this was the optimum. Possibly a noop.
        static_cast<derived_t const &>(*this).check_score(current_cell, cache);

        // Prepare horizontal and vertical score for next column.
        tmp += cache.gap_open_score;
        score_cell.up += cache.gap_extension_score;
        score_cell.w_left = score_cell.r_left + cache.gap_extension_score;

        score_cell.up = (score_cell.up < tmp) ? (trace_cell.up = trace_directions::up_open, tmp)
                                              : (trace_cell.up = trace_directions::up, score_cell.up);
        score_cell.w_left = (score_cell.w_left < tmp) ? (trace_cell.w_left = trace_directions::left_open, tmp)
                                                      : (trace_cell.w_left = trace_directions::left, score_cell.w_left);
    }

    /*!\brief Computes the score of the first cell within the band.
     * \tparam cell_t  The type of the current cell.
     * \tparam cache_t The type of the cache.
     * \param[in,out] current_cell    The current cell in the dynamic programming matrix.
     * \param[in,out] cache           The cache storing hot helper variables.
     * \param[in]     score           The score of comparing the respective letters of the first and the second sequence.
     */
    template <typename cell_t, typename cache_t>
    constexpr void compute_first_band_cell(cell_t && current_cell,
                                           cache_t & cache,
                                           score_t const score) const noexcept
    {
        // Compute the diagonal score and the compare with the previous horizontal value.
        auto [score_cell, trace_cell] = current_cell;
        score_cell.current = score_cell.diagonal + score;

        score_cell.current = (score_cell.current < score_cell.r_left)
                                ? (trace_cell.current = trace_cell.r_left, score_cell.r_left)
                                : (trace_cell.current = trace_directions::diagonal, score_cell.current);

        if constexpr (align_local_t::value)
        {
            score_cell.current = (score_cell.current < 0) ? (trace_cell.current = trace_directions::none, 0)
                                                          : score_cell.current;
        }
        // Check if this was the optimum. Possibly a noop.
        static_cast<derived_t const &>(*this).check_score(current_cell, cache);

        // At the top of the band we can not come from up but only diagonal or left, so the next vertical must be a
        // gap open.
        score_cell.up = score_cell.current + cache.gap_open_score;  // add gap open cost
        trace_cell.up = trace_directions::up_open;
    }

    /*!\brief Creates the cache used for affine gap computation.
     * \tparam    gap_scheme_t The type of the gap scheme.
     * \param[in] scheme       The configured gap scheme.
     * \returns The created cache.
     */
    template <typename gap_scheme_t>
    constexpr auto make_cache(gap_scheme_t && scheme) const noexcept
    {
        return alignment_algorithm_cache{static_cast<score_t>(scheme.get_gap_open_score() + scheme.get_gap_score()),
                                         static_cast<score_t>(scheme.get_gap_score())};
    }
};

} // namespace seqan3::detail
