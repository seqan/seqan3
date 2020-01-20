// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::affine_gap_init_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alignment/pairwise/detail/alignment_algorithm_state.hpp>

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

/*!\brief The CRTP-policy that implements the initialisation of the dynamic programming matrix with affine gaps.
 * \ingroup alignment_policy
 * \tparam alignment_algorithm_t The derived type (seqan3::detail::alignment_algorithm) to be augmented with this
 *                               CRTP-policy.
 * \tparam traits_type The traits type to determine the initialisation rules of the dynamic programming matrix.
 *                     Defaults to seqan3::detail::default_affine_init_traits.
 *
 * \remarks The template parameters of this CRTP-policy are selected in the
 *          seqan3::detail::alignment_configurator::select_gap_init_policy when selecting the alignment for the given
 *          configuration.
 */
template <typename alignment_algorithm_t, typename traits_type = default_affine_init_traits>
class affine_gap_init_policy
{
private:
    //!\brief Befriends the derived class to grant it access to the private members.
    friend alignment_algorithm_t;

    /*!\name Constructors, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr affine_gap_init_policy() noexcept = default;                                           //!< Defaulted
    constexpr affine_gap_init_policy(affine_gap_init_policy const &) noexcept = default;             //!< Defaulted
    constexpr affine_gap_init_policy(affine_gap_init_policy &&) noexcept = default;                  //!< Defaulted
    constexpr affine_gap_init_policy & operator=(affine_gap_init_policy const &) noexcept = default; //!< Defaulted
    constexpr affine_gap_init_policy & operator=(affine_gap_init_policy &&) noexcept = default;      //!< Defaulted
    ~affine_gap_init_policy() noexcept = default;                                                    //!< Defaulted
    //!\}

    /*!\brief Initialises the first cell of the dynamic programming matrix.
     * \tparam cell_t  The underlying cell type.
     * \tparam score_t The score type used inside of the alignment algorithm.
     *
     * \param[in,out] origin_cell The first cell of the dynamic programming matrix.
     * \param[in,out] state The state with gap information and the current alignment optimum.
     *
     * \details
     *
     * `cell_t` is the result type of dereferencing the zipped iterator over the respective alignment score matrix and
     * the alignment trace matrix used inside of the seqan3::detail::alignment_matrix_policy. The first parameter
     * stored in the zipped tuple is the seqan3::detail::alignment_score_matrix_proxy and the second value is the
     * seqan3::detail::alignment_trace_matrix_proxy.
     */
    template <typename cell_t, typename score_t>
    constexpr auto init_origin_cell(cell_t && origin_cell, alignment_algorithm_state<score_t> & state) const noexcept
    {
        // score_cell = seqan3::detail::alignment_score_matrix_proxy
        // trace_cell = seqan3::detail::alignment_trace_matrix_proxy
        auto & [score_cell, trace_cell] = origin_cell;

        // Initialise the first cell.
        score_cell.current = convert_to_simd_maybe<score_t>(0);
        trace_cell.current = convert_to_simd_maybe<score_t>(trace_directions::none);

        static_cast<alignment_algorithm_t const &>(*this).check_score_of_cell(origin_cell, state);

        // Initialise the vertical matrix cell according to the traits settings.
        if constexpr (traits_type::free_second_leading_t::value)
        {
            score_cell.up = convert_to_simd_maybe<score_t>(0);
            trace_cell.up = convert_to_simd_maybe<score_t>(trace_directions::none);
        }
        else // Initialise with gap_open score
        {
            score_cell.up = state.gap_open_score;
            trace_cell.up = convert_to_simd_maybe<score_t>(trace_directions::up_open);
        }

        // Initialise the horizontal matrix cell according to the traits settings.
        if constexpr (traits_type::free_first_leading_t::value)
        {
            score_cell.w_left = convert_to_simd_maybe<score_t>(0);
            trace_cell.w_left = convert_to_simd_maybe<score_t>(trace_directions::none);
        }
        else // Initialise with gap_open score
        {
            score_cell.w_left = state.gap_open_score;
            trace_cell.w_left = convert_to_simd_maybe<score_t>(trace_directions::left_open);
        }
    }

    /*!\brief Initialises a cell in the first column of the dynamic programming matrix.
     * \tparam cell_t  The underlying cell type.
     * \tparam score_t The score type used inside of the alignment algorithm.
     *
     * \param[in,out] column_cell A cell of the first row of the dynamic programming matrix.
     * \param[in,out] state The state with gap information and the current alignment optimum.
     *
     * \details
     *
     * `cell_t` is the result type of dereferencing the zipped iterator over the respective alignment score matrix and
     * the alignment trace matrix used inside of the seqan3::detail::alignment_matrix_policy. The first parameter
     * stored in the zipped tuple is the seqan3::detail::alignment_score_matrix_proxy and the second value is the
     * seqan3::detail::alignment_trace_matrix_proxy.
     */
    template <typename cell_t, typename score_t>
    constexpr auto init_column_cell(cell_t && column_cell, alignment_algorithm_state<score_t> & state) const noexcept
    {
        // score_cell = seqan3::detail::alignment_score_matrix_proxy
        // trace_cell = seqan3::detail::alignment_trace_matrix_proxy
        auto & [score_cell, trace_cell] = column_cell;

        score_cell.current = score_cell.up;
        trace_cell.current = trace_cell.up;

        static_cast<alignment_algorithm_t const &>(*this).check_score_of_cell(column_cell, state);

        // Initialise the vertical matrix cell according to the traits settings.
        if constexpr (traits_type::free_second_leading_t::value)
        {
            score_cell.up = convert_to_simd_maybe<score_t>(0);
        }
        else
        {
            score_cell.up += state.gap_extension_score;
            trace_cell.up = convert_to_simd_maybe<score_t>(trace_directions::up);
        }

        score_cell.w_left = score_cell.current + state.gap_open_score;
        trace_cell.w_left = convert_to_simd_maybe<score_t>(trace_directions::left_open);
    }

    /*!\brief Initialises a cell in the first row of the dynamic programming matrix.
     * \tparam cell_t  The underlying cell type.
     * \tparam score_t The score type used inside of the alignment algorithm.
     *
     * \param[in,out] row_cell A cell of the first row of the dynamic programming matrix.
     * \param[in,out] state The state with gap information and the current alignment optimum.
     *
     * \details
     *
     * `cell_t` is the result type of dereferencing the zipped iterator over the respective alignment score matrix and
     * the alignment trace matrix used inside of the seqan3::detail::alignment_matrix_policy. The first parameter
     * stored in the zipped tuple is the seqan3::detail::alignment_score_matrix_proxy and the second value is the
     * seqan3::detail::alignment_trace_matrix_proxy.
     */
    template <typename cell_t, typename score_t>
    constexpr auto init_row_cell(cell_t && row_cell, alignment_algorithm_state<score_t> & state) const noexcept
    {
        // score_cell = seqan3::detail::alignment_score_matrix_proxy
        // trace_cell = seqan3::detail::alignment_trace_matrix_proxy
        auto & [score_cell, trace_cell] = row_cell;

        score_cell.current = score_cell.r_left;
        trace_cell.current = trace_cell.r_left;

        static_cast<alignment_algorithm_t const &>(*this).check_score_of_cell(row_cell, state);

        score_cell.up = score_cell.current + state.gap_open_score;
        trace_cell.up = convert_to_simd_maybe<score_t>(trace_directions::up_open);

        if constexpr (traits_type::free_first_leading_t::value)
        {
            score_cell.w_left = convert_to_simd_maybe<score_t>(0);
            trace_cell.w_left = convert_to_simd_maybe<score_t>(trace_directions::none);
        }
        else
        {
            score_cell.w_left = score_cell.r_left + state.gap_extension_score;
            trace_cell.w_left = convert_to_simd_maybe<score_t>(trace_directions::left);
        }
    }

private:
    /*!\brief Converts the given value into a simd vector or just returns the value if alignment is not
     *        executed in vectorised mode.
     * \tparam score_t The type of the score; must model either seqan3::simd::simd_concept or
     *                 seqan3::arithmetic.
     * \tparam value_t The value type to convert; must model seqan3::arithmetic.
     *
     * \param[in] value The value to possibly convert.
     *
     * \returns A simd vector filled with the given value or the value itself if the alignment is not executed
     *          in vectorised mode.
     */
    template <typename score_t, typename value_t>
    constexpr auto convert_to_simd_maybe(value_t const value) const noexcept
    {
        if constexpr (simd_concept<score_t>)
        {
            using scalar_t = typename simd_traits<score_t>::scalar_type;
            return simd::fill<score_t>(static_cast<scalar_t>(value));
        }
        else
        {
            return value;
        }
    }

};

} // namespace seqan3::detail
