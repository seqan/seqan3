// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::affine_gap_init_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/core/algorithm/configuration.hpp>

namespace seqan3::detail
{

/*!\brief Implements the Initialisation of the dynamic programming matrix with affine gaps.
 * \ingroup alignment_policy
 * \tparam derived_t   The derived alignment algorithm.
 */
template <typename derived_t>
class affine_gap_init_policy
{
protected:
    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr affine_gap_init_policy() noexcept = default;
    constexpr affine_gap_init_policy(affine_gap_init_policy const &) noexcept = default;
    constexpr affine_gap_init_policy(affine_gap_init_policy &&) noexcept = default;
    constexpr affine_gap_init_policy & operator=(affine_gap_init_policy const &) noexcept = default;
    constexpr affine_gap_init_policy & operator=(affine_gap_init_policy &&) noexcept = default;
    ~affine_gap_init_policy() noexcept = default;
    //!}

    /*!\brief Initialises the origin of the dynamic programming matrix.
     * \tparam        cell_t       The underlying cell type.
     * \tparam        cache_t      The type of the cache.
     * \param[in,out] current_cell The current cell in the dynamic programming matrix.
     * \param[in,out] cache        The cache storing hot helper variables.
     */
    template <typename cell_t, typename cache_t>
    constexpr auto init_origin_cell(cell_t & current_cell, cache_t & cache) const noexcept
    {
        using std::get;

        auto & [main_score, hz_score] = current_cell;
        auto & [prev_cell, gap_open, gap_extend] = cache;
        auto & vt_score = get<1>(prev_cell);

        main_score = 0;
        hz_score = gap_open + gap_extend;
        vt_score = gap_open + gap_extend;
    }

    /*!\brief Initialises a cell in the first column of the dynamic programming matrix.
     * \tparam        cell_t       The underlying cell type.
     * \tparam        cache_t      The type of the cache.
     * \param[in,out] current_cell The current cell in the dynamic programming matrix.
     * \param[in,out] cache        The cache storing hot helper variables.
     */
    template <typename cell_t, typename cache_t>
    constexpr auto init_column_cell(cell_t & current_cell, cache_t & cache) const noexcept
    {
        using std::get;

        auto & [main_score, hz_score] = current_cell;
        auto & [prev_cell, gap_open, gap_extend] = cache;
        auto & vt_score = get<1>(prev_cell);

        main_score = vt_score;
        hz_score = main_score + gap_open + gap_extend;
        vt_score += gap_extend;
    }

    /*!\brief Initialises a cell in the first row of the dynamic programming matrix.
     * \tparam        cell_t       The underlying cell type.
     * \tparam        cache_t      The type of the cache.
     * \param[in,out] current_cell The current cell in the dynamic programming matrix.
     * \param[in,out] cache        The cache storing hot helper variables.
     */
    template <typename cell_t, typename cache_t>
    constexpr auto init_row_cell(cell_t & current_cell, cache_t & cache) const noexcept
    {
        auto & [main_score, hz_score] = current_cell;
        auto & [prev_cell, gap_open, gap_extend] = cache;
        auto & [prev_score, vt_score] = prev_cell;

        prev_score = main_score;
        main_score = hz_score;
        hz_score += gap_extend;
        vt_score += main_score + gap_open + gap_extend;
    }
};

} // namespace seqan3::detail
