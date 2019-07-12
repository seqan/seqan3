// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::affine_gap_banded_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <seqan3/alignment/pairwise/policy/affine_gap_policy.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// affine_gap_banded_policy
// ----------------------------------------------------------------------------

/*!\brief The hot kernel implementation using affine gaps.
 * \ingroup alignment_policy
 * \tparam derived_type   The derived alignment algorithm.
 * \tparam cell_type      The cell type of the dynamic programming matrix.
 * \tparam align_local_t  A bool constant that denotes whether local alignment is performed.
 */
template <typename derived_type, typename cell_type, typename align_local_t = std::false_type>
class affine_gap_banded_policy :
    public affine_gap_policy<affine_gap_banded_policy<derived_type, cell_type, align_local_t>, cell_type, align_local_t>
{
private:

    //!\brief The type of the base.
    using base_t = affine_gap_policy<affine_gap_banded_policy<derived_type, cell_type, align_local_t>,
                                     cell_type,
                                     align_local_t>;

    //!\brief Befriend the derived type.
    friend derived_type;

    //!\brief Import this typename into class scope.
    using score_t = typename base_t::score_t;

    /*!\name Constructors, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr affine_gap_banded_policy() noexcept = default;                                             //!< Defaulted
    constexpr affine_gap_banded_policy(affine_gap_banded_policy const &) noexcept = default;             //!< Defaulted
    constexpr affine_gap_banded_policy(affine_gap_banded_policy &&) noexcept = default;                  //!< Defaulted
    constexpr affine_gap_banded_policy & operator=(affine_gap_banded_policy const &) noexcept = default; //!< Defaulted
    constexpr affine_gap_banded_policy & operator=(affine_gap_banded_policy &&) noexcept = default;      //!< Defaulted
    ~affine_gap_banded_policy() noexcept = default;                                                      //!< Defaulted
    //!\}

    /*!\brief Computes the score for current cell.
     * \tparam        score_cell_type The type of the current cell.
     * \tparam        cache_t         The type of the cache.
     * \param[in,out] current_cell    The current cell in the dynamic programming matrix.
     * \param[in,out] cache           The cache storing hot helper variables.
     * \param[in]     score           The score of comparing the respective letters of the first and the second sequence.
     */
    template <typename score_cell_type, typename cache_t>
    constexpr void compute_cell(score_cell_type && current_cell,
                                cache_t & cache,
                                score_t const score) const noexcept
    {
        using std::get;
        // Unpack the cached variables.
        auto & [score_entry, coordinate, trace_value] = current_cell;
        auto & [current_entry, next_entry] = score_entry;
        auto & [main_score, hz_score, hz_trace] = current_entry;
        auto const & [prev_main_score, prev_hz_score, prev_hz_trace]= next_entry;
        auto & [prev_cell, gap_open, gap_extend, opt] = cache;
        auto & [tmp, vt_score, vt_trace] = prev_cell;

        (void) prev_main_score;

        //TODO For the local alignment we might be able to use GCC overlow builtin arithmetics, which
        // allows us to check if overflow/underflow would happen. Not sure, if this helps with the performance though.
        // (see https://gcc.gnu.org/onlinedocs/gcc/Integer-Overflow-Builtins.html)
        main_score += score;
        // Compute where the max comes from.
        if constexpr (decays_to_ignore_v<decltype(trace_value)>) // Don't compute a traceback
        {
            main_score = (main_score < vt_score) ? vt_score : main_score;
            main_score = (main_score < prev_hz_score) ? prev_hz_score : main_score;
            if constexpr (align_local_t::value)
                main_score = (main_score < 0) ? 0 : main_score;
        }
        else  // Compute any traceback
        {
            main_score = (main_score < vt_score) ? (trace_value = vt_trace, vt_score)
                                                 : (trace_value = trace_directions::diagonal | vt_trace, main_score);
            main_score = (main_score < prev_hz_score) ? (trace_value = prev_hz_trace, prev_hz_score)
                                                      : (trace_value |= prev_hz_trace, main_score);
            if constexpr (align_local_t::value)
                main_score = (main_score < 0) ? (trace_value = trace_directions::none, 0) : main_score;
        }

        // Check if this was the optimum. Possibly a noop.
        static_cast<derived_type const &>(*this).check_score(
            alignment_optimum{main_score, static_cast<alignment_coordinate>(coordinate)}, opt);

        tmp = main_score + gap_open;
        vt_score += gap_extend;
        hz_score = prev_hz_score + gap_extend;

        if constexpr (decays_to_ignore_v<decltype(trace_value)>)
        {
            vt_score = (vt_score < tmp) ? tmp : vt_score;
            hz_score = (hz_score < tmp) ? tmp : hz_score;
        }
        else
        {
            vt_score = (vt_score < tmp) ? (vt_trace = trace_directions::up_open, tmp)
                                        : (vt_trace = trace_directions::up, vt_score);
            hz_score = (hz_score < tmp) ? (hz_trace = trace_directions::left_open, tmp)
                                        : (hz_trace = trace_directions::left, hz_score);
        }
    }

    /*!\brief Computes the score of the first cell within the band.
     * \tparam        score_cell_type The type of the current cell.
     * \tparam        cache_t         The type of the cache.
     * \param[in,out] current_cell    The current cell in the dynamic programming matrix.
     * \param[in,out] cache           The cache storing hot helper variables.
     * \param[in]     score           The score of comparing the respective letters of the first and the second sequence.
     */
    template <typename score_cell_type, typename cache_t>
    constexpr void compute_first_band_cell(score_cell_type && current_cell,
                                           cache_t & cache,
                                           score_t const score) const noexcept
    {
        using std::get;
        // Unpack the cached variables.
        auto & [score_entry, coordinate, trace_value] = current_cell;
        auto & [current_entry, next_entry] = score_entry;
        auto & main_score = get<0>(current_entry);
        auto const & prev_hz_score = get<1>(next_entry);
        auto const & prev_hz_trace = get<2>(next_entry);
        auto & [tmp, vt_score, vt_trace] = get<0>(cache);

        (void) tmp;

        // Compute the diagonal score and the compare with the previous horizontal value.
        main_score += score;

        if constexpr (decays_to_ignore_v<decltype(trace_value)>) // Don't compute a traceback
        {
            main_score = (main_score < prev_hz_score) ? prev_hz_score : main_score;
            if constexpr (align_local_t::value)
                main_score = (main_score < 0) ? 0 : main_score;
        }
        else
        {
            main_score = (main_score < prev_hz_score) ? (trace_value = prev_hz_trace, prev_hz_score)
                                   : (trace_value = trace_directions::diagonal, main_score);
            if constexpr (align_local_t::value)
                main_score = (main_score < 0) ? (trace_value = trace_directions::none, 0) : main_score;
        }

        // Check if this was the optimum. Possibly a noop.
        static_cast<derived_type const &>(*this).check_score(
                alignment_optimum{main_score, static_cast<alignment_coordinate>(coordinate)}, get<3>(cache));
        // At the top of the band we can not come from up but only diagonal or left, so the next vertical must be a
        // gap open.
        vt_score = main_score + get<1>(cache);  // add gap open cost

        // Store vertical open as it is the only way to get into the first cell of the band from a vertical direction.
        if constexpr (!decays_to_ignore_v<decltype(trace_value)>) // Traceback requested.
            vt_trace = trace_directions::up_open;

        // Don't compute the horizontal value, since it will never be used.
    }

    using base_t::make_cache;
};

} // namespace seqan3::detail
