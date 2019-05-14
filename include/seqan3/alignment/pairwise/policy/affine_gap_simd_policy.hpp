// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::affine_gap_simd_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <limits>
#include <tuple>

#include <seqan3/alignment/matrix/alignment_optimum.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// affine_gap_simd_policy
// ----------------------------------------------------------------------------

/*!\brief The hot kernel implementation using affine gaps.
 * \ingroup alignment_policy
 * \tparam derived_t   The derived alignment algorithm.
 * \tparam cell_t      The cell type of the dynamic programming matrix.
 */
template <typename derived_t, typename cell_t>
class affine_gap_simd_policy
{
private:

    //!\brief Befriends the derived class to grant it access to the private members.
    friend derived_t;

    /*!\name Member types
     * \{
     */
    //!\brief The underlying score type.
    using score_t = std::tuple_element_t<0, cell_t>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr affine_gap_simd_policy() noexcept = default; //!<Defaulted.
    constexpr affine_gap_simd_policy(affine_gap_simd_policy const &) noexcept = default; //!<Defaulted.
    constexpr affine_gap_simd_policy(affine_gap_simd_policy &&) noexcept = default; //!<Defaulted.
    constexpr affine_gap_simd_policy & operator=(affine_gap_simd_policy const &) noexcept = default; //!<Defaulted.
    constexpr affine_gap_simd_policy & operator=(affine_gap_simd_policy &&) noexcept = default; //!<Defaulted.
    ~affine_gap_simd_policy() noexcept = default; //!<Defaulted.
    //!\}

    /*!\brief Computes the score of the current cell.
     * \tparam  matrix_entry_type The type of the current cell.
     * \tparam  cache_type        The type of the cache.
     * \param[in,out] current_cell The current cell in the dynamic programming matrix.
     * \param[in,out] cache        The cache storing hot helper variables.
     * \param[in]     score        The score of comparing the respective letters of the first and the second sequence.
     */
    template <typename matrix_entry_type, typename cache_t>
    constexpr void compute_cell(matrix_entry_type && current_cell,
                                cache_t & cache,
                                score_t const score) const noexcept
    {
        using std::get;

        // Unpack the cached variables.
        auto & [score_entry, coordinate, trace_value] = current_cell;
        auto & [main_score, hz_score, hz_trace]       = score_entry;
        auto & [prev_cell, gap_open, gap_extend, opt] = cache;
        auto & [prev_score, vt_score, vt_trace]       = prev_cell;

        // Precompute the diagonal score.
        score_t tmp = prev_score + score;

        // Compute where the max comes from.
        if constexpr (decays_to_ignore_v<decltype(trace_value)>) // Don't compute a traceback
        {
            tmp = (tmp < vt_score) ? vt_score : tmp;
            tmp = (tmp < hz_score) ? hz_score : tmp;
        }
        else  // Compute any traceback
        {
            // tmp = (tmp < vt_score) ? (trace_value = vt_trace, vt_score)
            //                        : (trace_value = trace_directions::diagonal | vt_trace, tmp);
            // tmp = (tmp < hz_score) ? (trace_value = hz_trace, hz_score)
            //                        : (trace_value |= hz_trace, tmp);
        }
        // Cache the current main score for the next diagonal computation and update the current score.
        prev_score = main_score;
        main_score = tmp;
        // Check if this was the optimum. Possibly a noop.
        static_cast<derived_t const &>(*this).check_score(
                alignment_optimum{tmp, static_cast<alignment_coordinate>(coordinate)}, opt);

        // Prepare horizontal and vertical score for next column.
        tmp += gap_open;
        vt_score += gap_extend;
        hz_score += gap_extend;

        if constexpr (decays_to_ignore_v<decltype(trace_value)>)
        {
            vt_score = (vt_score < tmp) ? tmp : vt_score;
            hz_score = (hz_score < tmp) ? tmp : hz_score;
        }
        else
        {
            // vt_score = (vt_score < tmp) ? (vt_trace = trace_directions::up_open, tmp)
            //                             : (vt_trace = trace_directions::up, vt_score);
            // hz_score = (hz_score < tmp) ? (hz_trace = trace_directions::left_open, tmp)
            //                             : (hz_trace = trace_directions::left, hz_score);
        }
    }

    /*!\brief Creates the cache used for affine gap computation.
     * \tparam    gap_scheme_t The type of the gap scheme.
     * \param[in] scheme       The configured gap scheme.
     * \returns The created cache.
     */
    template <typename gap_scheme_t>
    constexpr auto make_cache(gap_scheme_t && scheme) const noexcept
    {
        alignment_optimum opt{std::numeric_limits<score_t>::lowest(),
                              alignment_coordinate{column_index_type{0u}, row_index_type{0u}}};
        return std::tuple{cell_t{},
                          static_cast<score_t>(scheme.get_gap_open_score() + scheme.get_gap_score()),
                          static_cast<score_t>(scheme.get_gap_score()),
                          std::move(opt)};
    }
};

} // namespace seqan3::detail
