// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::affine_gap_banded_init_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/alignment/configuration/align_config_aligned_ends.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_init_policy.hpp>
#include <seqan3/core/algorithm/configuration.hpp>

namespace seqan3::detail
{

/*!\brief Implements the Initialisation of the dynamic programming matrix with affine gaps.
 * \ingroup alignment_policy
 * \tparam derived_t   The derived alignment algorithm.
 * \tparam traits_type The traits type to determine the initialisation rules of the dynamic programming matrix.
 *                     Defaults to seqan3::detail::default_affine_init_traits.
 */
template <typename derived_t, typename traits_type = default_affine_init_traits>
class affine_gap_banded_init_policy :
    public affine_gap_init_policy<affine_gap_banded_init_policy<derived_t, traits_type>, traits_type>
{
private:

    //!\brief The base type.
    using base_t = affine_gap_init_policy<affine_gap_banded_init_policy<derived_t, traits_type>, traits_type>;

    //!\brief Befriends the derived type.
    friend derived_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr affine_gap_banded_init_policy() noexcept = default;                                      //!< Defaulted
    constexpr affine_gap_banded_init_policy(affine_gap_banded_init_policy const &) noexcept = default; //!< Defaulted
    constexpr affine_gap_banded_init_policy(affine_gap_banded_init_policy &&) noexcept = default;      //!< Defaulted
    //!\brief Defaulted
    constexpr affine_gap_banded_init_policy & operator=(affine_gap_banded_init_policy const &) noexcept = default;
    //!\brief Defaulted
    constexpr affine_gap_banded_init_policy & operator=(affine_gap_banded_init_policy &&) noexcept = default;
    ~affine_gap_banded_init_policy() noexcept = default;                                               //!< Defaulted
    //!\}

    /*!\brief Initialises the origin of the dynamic programming matrix.
     * \tparam        cell_t       The underlying cell type.
     * \tparam        cache_t      The type of the cache.
     * \param[in,out] current_cell The current cell in the dynamic programming matrix.
     * \param[in,out] cache        The cache storing hot helper variables.
     */
    template <typename cell_t, typename cache_t>
    constexpr auto init_origin_cell(cell_t && current_cell, cache_t & cache) const noexcept
    {
        using std::get;

        // Call get twice since the banded score column is a zipped range and we need to access the main column.
        auto & [main_score, hz_score, hz_trace] = get<0>(get<0>(current_cell));
        auto & trace_value = get<2>(current_cell);
        auto & vt_score = get<1>(get<0>(cache));
        auto & vt_trace = get<2>(get<0>(cache));

        main_score = 0;

        // Initialise the vertical matrix cell according to the traits settings.
        if constexpr (traits_type::free_second_leading_t::value)
        {
            vt_score = 0;
            vt_trace = trace_directions::none;
        }
        else
        {
            vt_score = get<1>(cache);
            vt_trace = trace_directions::up_open;
        }
        // Initialise the horizontal matrix cell according to the traits settings.
        if constexpr (traits_type::free_first_leading_t::value)
        {
            hz_score = 0;
            hz_trace = trace_directions::none;
        }
        else
        {
            hz_score = get<1>(cache);
            hz_trace = trace_directions::left_open;
        }
        trace_value = trace_directions::none;
    }

    /*!\brief Initialises a cell in the first column of the dynamic programming matrix.
     * \tparam        cell_t       The underlying cell type.
     * \tparam        cache_t      The type of the cache.
     * \param[in,out] current_cell The current cell in the dynamic programming matrix.
     * \param[in,out] cache        The cache storing hot helper variables.
     */
    template <typename cell_t, typename cache_t>
    constexpr auto init_column_cell(cell_t && current_cell, cache_t & cache) const noexcept
    {
        using std::get;

        // Call get twice since the banded score column is a zipped range and we need to access the main column.
        auto & [main_score, hz_score, hz_trace] = get<0>(get<0>(current_cell));
        auto & trace_value = get<2>(current_cell);
        auto & vt_score = get<1>(get<0>(cache));
        auto & vt_trace = get<2>(get<0>(cache));

        main_score = vt_score;
        trace_value = vt_trace;

        // Initialise the vertical matrix cell according to the traits settings.
        if constexpr (traits_type::free_second_leading_t::value)
        {
            vt_score = 0;
            vt_trace = trace_directions::none;
        }
        else
        { // previous vertical + gap extension
            vt_score += get<2>(cache);
            vt_trace = trace_directions::up;
        }
        hz_score = main_score + get<1>(cache); // gap_opening cost
        hz_trace = trace_directions::left_open;
    }

    /*!\brief Initialises a cell in the first row of the current band.
     * \tparam        cell_t       The underlying cell type.
     * \tparam        cache_t      The type of the cache.
     * \param[in,out] current_cell The current cell in the dynamic programming matrix.
     * \param[in,out] cache        The cache storing hot helper variables.
     */
    template <typename cell_t, typename cache_t>
    constexpr auto init_row_cell(cell_t && current_cell, cache_t & cache) const noexcept
    {
        auto & [current_entry, next_entry] = get<0>(current_cell); // Split into current entry and next entry.
        auto & trace_value = get<2>(current_cell);  // the trace value to store.
        auto & main_score = get<0>(current_entry);  // current_entry stores current score to be updated.
        auto & hz_trace = get<2>(current_entry); // store trace for the next horizontal computation.
        auto const & prev_hz_score = get<1>(next_entry);  // get the previous horizontal score.
        auto const & prev_hz_trace = get<2>(next_entry);  // get the previous horizontal trace.
        auto & vt_score = get<1>(get<0>(cache));
        auto & vt_trace = get<2>(get<0>(cache));

        main_score = prev_hz_score;
        trace_value = prev_hz_trace;
        vt_score += main_score + get<1>(cache); // gap opening cost
        vt_trace = trace_directions::up_open;

        // Initialise the horizontal matrix cell according to the traits settings.
        if constexpr (traits_type::free_first_leading_t::value)
        {
            get<1>(current_entry) = 0;
            hz_trace = trace_directions::none;
        }
        else
        {
            get<1>(current_entry) = prev_hz_score + get<2>(cache);
            hz_trace = trace_directions::left;
        }
    }

    /*!\brief Balances the total score based on the band parameters and the alignment configuration.
     * \tparam        optimum_type     The type of the optimum to update.
     * \tparam        band_type        The type of the band.
     * \tparam        gap_scheme_type  The type of the gap_scheme.
     * \param[in,out] total  The total score to update.
     * \param[in]     band   The band settings.
     * \param[in]     scheme The gap scheme to compute the additional gap score.
     *
     * \details
     *
     * Depending on the band position and the alignment configuration updates the total score
     * of the alignment. It adds the score for initialising the matrix with a gap until the begin of the band.
     */
    template <typename optimum_type, typename band_type, typename gap_scheme_type>
    constexpr void balance_leading_gaps(optimum_type & total,
                                        band_type const & band,
                                        gap_scheme_type const & scheme) const noexcept
    {
        if constexpr (!traits_type::free_second_leading_t::value)
        {  // Band starts inside of second sequence.
            if (0 > band.upper_bound)
                total.score += scheme.score(std::abs(band.upper_bound));
        }
        if constexpr (!traits_type::free_first_leading_t::value)
        { // Band starts inside of first sequence.
            if (band.lower_bound > 0)
                total.score += scheme.score(band.lower_bound);
        }
    }
};

} // namespace seqan3::detail
