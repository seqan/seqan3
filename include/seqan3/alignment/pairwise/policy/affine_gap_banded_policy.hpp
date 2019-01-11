// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::affine_gap_banded_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <seqan3/alignment/pairwise/policy/affine_gap_policy.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/range.hpp>
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
 */
template <typename derived_type, typename cell_type>
class affine_gap_banded_policy : public affine_gap_policy<affine_gap_banded_policy<derived_type, cell_type>, cell_type>
{
private:

    //!\brief The type of the base.
    using base_t = affine_gap_policy<affine_gap_banded_policy<derived_type, cell_type>, cell_type>;

    //!\brief Befriend the derived type.
    friend derived_type;

    //!\brief Import this typename into class scope.
    using score_t = typename base_t::score_t;

    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr affine_gap_banded_policy() noexcept = default;
    constexpr affine_gap_banded_policy(affine_gap_banded_policy const &) noexcept = default;
    constexpr affine_gap_banded_policy(affine_gap_banded_policy &&) noexcept = default;
    constexpr affine_gap_banded_policy & operator=(affine_gap_banded_policy const &) noexcept = default;
    constexpr affine_gap_banded_policy & operator=(affine_gap_banded_policy &&) noexcept = default;
    ~affine_gap_banded_policy() noexcept = default;
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
        auto & [current_entry, next_entry] = get<0>(current_cell);
        auto & [main_score, hz_score] = current_entry;
        auto const & prev_hz_score = get<1>(next_entry);
        auto & [prev_cell, gap_open, gap_extend, opt] = cache;
        auto & [tmp, vt_score] = prev_cell;

        //TODO For the local alignment we might be able to use GCC overlow builtin arithmetics, which
        // allows us to check if overflow/underflow would happen. Not sure, if this helps with the performance though.
        // (see https://gcc.gnu.org/onlinedocs/gcc/Integer-Overflow-Builtins.html)
        main_score += score;
        main_score = std::max(main_score, vt_score);
        main_score = std::max(main_score, prev_hz_score);
        // Check if this was the optimum. Possibly a noop.
        static_cast<derived_type const &>(*this).check_score(main_score, opt);

        tmp = main_score + gap_open;
        vt_score = std::max(vt_score + gap_extend, tmp);
        hz_score = std::max(prev_hz_score + gap_extend, tmp);
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
        auto & [current_entry, next_entry] = get<0>(current_cell);
        auto & main_score = get<0>(current_entry);
        auto const & prev_hz_score = get<1>(next_entry);
        auto & vt_score = get<1>(get<0>(cache));

        // Compute the diagonal score and the compare with the previous horizontal value.
        main_score += score;
        main_score = std::max(main_score, prev_hz_score);
        // Check if this was the optimum. Possibly a noop.
        static_cast<derived_type const &>(*this).check_score(main_score, get<3>(cache));
        // At the top of the band we can not come from up but only diagonal or left, so the next vertical must be a
        // gap open.
        vt_score = main_score + get<1>(cache);  // add gap open cost
        // Don't compute the horizontal value, as it won't be used anyway.
    }

    using base_t::make_cache;
};

} // namespace seqan3::detail
