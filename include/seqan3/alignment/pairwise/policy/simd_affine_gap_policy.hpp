// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::simd_affine_gap_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <limits>
#include <tuple>

#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/matrix/alignment_optimum.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alignment/pairwise/detail/alignment_algorithm_state.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/simd_traits.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// simd_affine_gap_policy
// ----------------------------------------------------------------------------

/*!\brief The CRTP-policy that computes a batch of cells in the alignment matrix using simd instructions.
 * \ingroup alignment_policy
 * \tparam alignment_algorithm_t The derived type (seqan3::detail::alignment_algorithm) to be augmented with this
 *                               CRTP-policy.
 * \tparam score_t The score type of the dynamic programming matrix; must model seqan3::simd::simd_concept.
 * \tparam align_local_t A std::bool_constant to switch between local and global alignment.
 *
 * \details
 *
 * This CRTP-policy implements the recursion for the alignment algorithm with affine gaps using an inter-sequence
 * vectorisation scheme. See `Rahn, R, et al. Generic accelerated sequence alignment in SeqAn using vectorization
 * and multi-threading. Bioinformatics 34.20 (2018): 3437-3445.` for more information.
 *
 * \remarks The template parameters of this CRTP-policy are selected in the
 *          seqan3::detail::alignment_configurator::select_gap_policy when selecting the alignment for the given
 *          configuration.
 */
template <typename alignment_algorithm_t, simd_concept score_t, typename align_local_t = std::false_type>
class simd_affine_gap_policy
{
private:
    //!\brief Befriends the derived class to grant it access to the private members.
    friend alignment_algorithm_t;

    //!\brief The type of state of the alignment algorithm for affine gaps.
    using alignment_state_t = alignment_algorithm_state<score_t>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr simd_affine_gap_policy() noexcept = default; //!< Defaulted.
    constexpr simd_affine_gap_policy(simd_affine_gap_policy const &) noexcept = default; //!< Defaulted.
    constexpr simd_affine_gap_policy(simd_affine_gap_policy &&) noexcept = default; //!< Defaulted.
    constexpr simd_affine_gap_policy & operator=(simd_affine_gap_policy const &) noexcept = default; //!< Defaulted.
    constexpr simd_affine_gap_policy & operator=(simd_affine_gap_policy &&) noexcept = default; //!< Defaulted.
    ~simd_affine_gap_policy() noexcept = default; //!< Defaulted.
    //!\}

    /*!\brief Computes the score of the current simd cell.
     * \tparam cell_t The type of the current cell [for detailed information on the type see below].
     * \param[in,out] current_cell The current cell in the dynamic programming matrix.
     * \param[in,out] state        The state storing hot helper variables.
     * \param[in]     score        The score of comparing the respective letters of the first and the second sequence.
     *
     * \details
     *
     * `cell_t` is the result type of dereferencing the zipped iterator over the respective alignment score matrix and
     * the alignment trace matrix used inside of the seqan3::detail::alignment_matrix_policy. The first parameter
     * stored in the zipped tuple is the seqan3::detail::alignment_score_matrix_proxy and the second value is the
     * seqan3::detail::alignment_trace_matrix_proxy.
     *
     * In order to compute the maximum for two simd vectors, gcc implements the ternary operator such that
     * `std::max(a, b)` can be implemented as `(a > b) ? a : b`, where `(a > b)` returns a mask vector. This implements
     * the compare-and-blend approach for simd vector types.
     */
    template <typename cell_t>
    constexpr void compute_cell(cell_t && current_cell,
                                alignment_algorithm_state<score_t> & state,
                                score_t const score) const noexcept
    {
        // score_cell = seqan3::detail::alignment_score_matrix_proxy
        // trace_cell = seqan3::detail::alignment_trace_matrix_proxy
        auto & [score_cell, trace_cell] = current_cell;
        constexpr bool with_trace = !decays_to_ignore_v<std::remove_reference_t<decltype(trace_cell.current)>>;
        // Precompute the diagonal score.
        score_t tmp = score_cell.diagonal + score;

        if constexpr (with_trace)
        {
            auto mask = tmp < score_cell.up;
            tmp = (mask) ? score_cell.up : tmp;
            trace_cell.current = (mask) ? trace_cell.up : convert_to_simd(trace_directions::diagonal) | trace_cell.up;

            mask = tmp < score_cell.r_left;
            tmp = (mask) ? score_cell.r_left : tmp;
            trace_cell.current = (mask) ? trace_cell.r_left : trace_cell.current | trace_cell.r_left;
        }
        else
        {
            tmp = (tmp < score_cell.up) ? score_cell.up : tmp;
            tmp = (tmp < score_cell.r_left) ? score_cell.r_left : tmp;
        }

        if constexpr (align_local_t::value)
        {
            tmp = (tmp < simd::fill<score_t>(0))
            /*then*/ ? (trace_cell.current = convert_to_simd(trace_directions::none), simd::fill<score_t>(0))
            /*else*/ : tmp;
        }

        // Store the current max score.
        score_cell.current = tmp;
        // Check if this was the optimum. Possibly a noop.
        static_cast<alignment_algorithm_t const &>(*this).check_score_of_cell(current_cell, state);

        // Prepare horizontal and vertical score for next column.
        tmp += state.gap_open_score;
        score_cell.up += state.gap_extension_score;
        score_cell.w_left = score_cell.r_left + state.gap_extension_score;

        auto mask = score_cell.up < tmp;
        score_cell.up = (mask) ? tmp : score_cell.up;
        trace_cell.up = (mask) ? convert_to_simd(trace_directions::up_open) : convert_to_simd(trace_directions::up);
        mask = score_cell.w_left < tmp;
        score_cell.w_left = (mask) ? tmp : score_cell.w_left;
        trace_cell.w_left = (mask) ? convert_to_simd(trace_directions::left_open)
                                   : convert_to_simd(trace_directions::left);
    }

    /*!\brief Initialise the alignment state for affine gap computation.
     * \tparam alignment_configuration_t The type of alignment configuration.
     * \param[in] config The alignment configuration.
     * \returns The initialised seqan3::detail::alignment_algorithm_state.
     *
     * \details
     *
     * Gets the stored gap scheme from the configuration or uses a default gap scheme if not set and initialises the
     * alignment algorithm state with the respective gap extension and gap open costs. If the gap scheme was not
     * specified by the user the following defaults are used:
     *  * `-1` for the gap extension score, and
     *  * `-10` for the gap open score.
     */
    template <typename alignment_configuration_t>
    constexpr void initialise_alignment_state(alignment_configuration_t const & config) noexcept
    {
        using scalar_t = typename simd_traits<score_t>::scalar_type;
        auto scheme = config.template value_or<align_cfg::gap>(gap_scheme{gap_score{-1}, gap_open_score{-10}});

        alignment_state.gap_extension_score = simd::fill<score_t>(static_cast<scalar_t>(scheme.get_gap_score()));
        alignment_state.gap_open_score = simd::fill<score_t>(static_cast<scalar_t>(scheme.get_gap_score() +
                                                                                   scheme.get_gap_open_score()));
    }

    /*!\brief Converts a trace direction into a simd vector.
     * \param[in] direction The trace direction to convert to a simd vector.
     */
    constexpr score_t convert_to_simd(trace_directions const direction) const noexcept
    {
        using scalar_t = typename simd_traits<score_t>::scalar_type;

        return simd::fill<score_t>(static_cast<scalar_t>(direction));
    }

    alignment_state_t alignment_state{}; //!< The internal alignment state tracking the current alignment optimum.
};

} // namespace seqan3::detail
