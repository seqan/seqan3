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

#include <seqan3/alignment/matrix/trace_directions.hpp>

namespace seqan3::detail
{

/*!\brief The default traits class for seqan3::detail::affine_gap_init_policy.
 * \ingroup alignment_policy
 *
 * \details
 *
 * Enables the behaviour of a global alignment where both sides of the dynamic programming matrix are
 * initialised with growing gap penalties.
 */
struct default_affine_init_traits
{
    //!\brief A std::bool_constant to enable/disable free end gaps for the leading gaps of the first sequence.
    using free_first_leading_t  = std::false_type;
    //!\brief A std::bool_constant to enable/disable free end gaps for the leading gaps of the second sequence.
    using free_second_leading_t = std::false_type;
};

/*!\brief Implements the initialisation of the dynamic programming matrix with affine gaps.
 * \ingroup alignment_policy
 * \tparam derived_t   The derived alignment algorithm.
 * \tparam traits_type The traits type to determine the initialisation rules of the dynamic programming matrix.
 *                     Defaults to seqan3::detail::default_affine_init_traits.
 */
template <typename derived_t, typename traits_type = default_affine_init_traits>
class affine_gap_init_policy
{
private:

    //!\brief Befriends the derived class to grant it access to the private members.
    friend derived_t;
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

        auto & [main_score, hz_score, hz_trace] = get<0>(current_cell);
        auto & prev_cell = get<0>(cache);
        auto & vt_score = get<1>(prev_cell);

        main_score = 0;
        get<2>(current_cell) = trace_directions::none; // store the trace direction

        // Initialise the vertical matrix cell according to the traits settings.
        if constexpr (traits_type::free_second_leading_t::value)
        {
            vt_score = 0;
            get<2>(prev_cell) = trace_directions::none;  // cache vertical trace
        }
        else
        {
            vt_score = get<1>(cache);
            get<2>(prev_cell) = trace_directions::up_open; // cache vertical trace
        }

        // Initialise the horizontal matrix cell according to the traits settings.
        if constexpr (traits_type::free_first_leading_t::value)
        {
            hz_score = 0;
            hz_trace = trace_directions::none; // cache horizontal trace
        }
        else
        {
            hz_score = get<1>(cache);
            hz_trace = trace_directions::left_open; // cache horizontal trace
        }
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

        auto & [main_score, hz_score, hz_trace] = get<0>(current_cell);
        auto & prev_cell = get<0>(cache);
        auto & vt_score = get<1>(prev_cell);

        main_score = vt_score;
        get<2>(current_cell) = get<2>(prev_cell); // store the trace direction

        // Initialise the vertical matrix cell according to the traits settings.
        if constexpr (traits_type::free_second_leading_t::value)
        {
            vt_score = 0;
            get<2>(prev_cell) = trace_directions::none; // cache vertical trace
        }
        else
        {
            vt_score += get<2>(cache);
            get<2>(prev_cell) = trace_directions::up; // cache vertical trace
        }

        hz_score = main_score + get<1>(cache);
        hz_trace = trace_directions::none;  // cache horizontal trace
    }

    /*!\brief Initialises a cell in the first row of the dynamic programming matrix.
     * \tparam        cell_t       The underlying cell type.
     * \tparam        cache_t      The type of the cache.
     * \param[in,out] current_cell The current cell in the dynamic programming matrix.
     * \param[in,out] cache        The cache storing hot helper variables.
     */
    template <typename cell_t, typename cache_t>
    constexpr auto init_row_cell(cell_t && current_cell, cache_t & cache) const noexcept
    {
        using std::get;

        auto & [main_score, hz_score, hz_trace] = get<0>(current_cell);
        auto & [prev_score, vt_score, vt_trace] = get<0>(cache);

        prev_score = main_score;
        main_score = hz_score;
        get<2>(current_cell) = hz_trace; // store the trace direction

        vt_score += main_score + get<1>(cache);
        vt_trace = trace_directions::none; // cache vertical trace
        // Initialise the horizontal matrix cell according to the traits settings.
        if constexpr (traits_type::free_first_leading_t::value)
        {
            hz_score = 0;
            hz_trace = trace_directions::none;
        }
        else
        {
            hz_score += get<2>(cache);
            hz_trace = trace_directions::left;
        }
    }
};

} // namespace seqan3::detail
