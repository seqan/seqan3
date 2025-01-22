// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::affine_gap_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/alignment/matrix/detail/trace_directions.hpp>
#include <seqan3/alignment/pairwise/detail/alignment_algorithm_state.hpp>
#include <seqan3/utility/concept.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// affine_gap_policy
// ----------------------------------------------------------------------------

/*!\brief The CRTP-policy that computes a single cell in the alignment matrix.
 * \ingroup alignment_pairwise_policy
 * \tparam alignment_algorithm_t The derived type (seqan3::detail::alignment_algorithm) to be augmented with this
 *                               CRTP-policy.
 * \tparam score_t The score type of the dynamic programming matrix; must model seqan3::arithmetic.
 * \tparam align_local_t A std::bool_constant to switch between local and global alignment.
 *
 * \details
 *
 * This CRTP-policy implements the recursion for the alignment algorithm with affine gaps. See
 * `Osamu Gotoh, Optimal sequence alignment allowing for long gaps, Bulletin of Mathematical Biology, Volume 52,
 *  Issue 3, 1990, Pages 359-373, ISSN 0092-8240, https://doi.org/10.1007/BF02458577`
 * for more information.
 *
 * \remarks The template parameters of this CRTP-policy are selected in the
 *          seqan3::detail::alignment_configurator::select_gap_policy when selecting the alignment for the given
 *          configuration.
 */
template <typename alignment_algorithm_t, typename score_t, typename align_local_t = std::false_type>
class affine_gap_policy
{
private:
    //!\brief Befriends the derived class to grant it access to the private members.
    friend alignment_algorithm_t;

    //!\brief The state of the alignment algorithm for affine gaps.
    using alignment_state_t = alignment_algorithm_state<score_t>;

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

    //!\brief Initialise the policy.
    template <typename configuration_t>
    affine_gap_policy(configuration_t const & /*config*/)
    {}
    //!\}

    /*!\brief Computes the score of the current cell.
     * \tparam cell_t The type of the current cell [for detailed information on the type see below].
     * \param[in,out] current_cell The current cell in the dynamic programming matrix.
     * \param[in,out] cache        The cache storing hot helper variables.
     * \param[in]     score        The score of comparing the respective letters of the first and the second sequence.
     *
     * \details
     *
     * `cell_t` is the result type of dereferencing the zipped iterator over the respective alignment score matrix and
     * the alignment trace matrix used inside of the seqan3::detail::alignment_matrix_policy. The first parameter
     * stored in the zipped tuple is the seqan3::detail::alignment_score_matrix_proxy and the second value is the
     * seqan3::detail::alignment_trace_matrix_proxy.
     */
    template <typename cell_t>
    constexpr void
    compute_cell(cell_t && current_cell, alignment_algorithm_state<score_t> & cache, score_t const score) const noexcept
    {
        // score_cell = seqan3::detail::alignment_score_matrix_proxy
        // trace_cell = seqan3::detail::alignment_trace_matrix_proxy
        auto & [score_cell, trace_cell] = current_cell;
        constexpr bool with_trace = !decays_to_ignore_v<std::remove_reference_t<decltype(trace_cell.current)>>;
        // Precompute the diagonal score.
        score_t tmp = score_cell.diagonal + score;

        if constexpr (with_trace)
        {
            tmp = (tmp < score_cell.up) ? (trace_cell.current = trace_cell.up, score_cell.up)
                                        : (trace_cell.current = trace_directions::diagonal | trace_cell.up, tmp);
            tmp = (tmp < score_cell.r_left)
                    ? (trace_cell.current = trace_cell.r_left | (trace_cell.current & trace_directions::carry_up_open),
                       score_cell.r_left)
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
        static_cast<alignment_algorithm_t const &>(*this).check_score_of_cell(current_cell, cache);

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
     * \tparam cell_t The type of the current cell [for detailed information on the type see below].
     * \param[in,out] current_cell    The current cell in the dynamic programming matrix.
     * \param[in,out] cache           The cache storing hot helper variables.
     * \param[in]     score           The score of comparing the respective letters of the first and the second sequence.
     *
     * \details
     *
     * `cell_t` is the result type of dereferencing the zipped iterator over the respective alignment score matrix and
     * the alignment trace matrix used inside of the seqan3::detail::alignment_matrix_policy. The first parameter
     * stored in the zipped tuple is the seqan3::detail::alignment_score_matrix_proxy and the second value is the
     * seqan3::detail::alignment_trace_matrix_proxy.
     */
    template <typename cell_t>
    constexpr void compute_first_band_cell(cell_t && current_cell,
                                           alignment_algorithm_state<score_t> & cache,
                                           score_t const score) const noexcept
    {
        // Compute the diagonal score and the compare with the previous horizontal value.
        // score_cell = seqan3::detail::alignment_score_matrix_proxy
        // trace_cell = seqan3::detail::alignment_trace_matrix_proxy
        auto & [score_cell, trace_cell] = current_cell;
        score_cell.current = score_cell.diagonal + score;

        score_cell.current = (score_cell.current < score_cell.r_left)
                               ? (trace_cell.current = trace_cell.r_left, score_cell.r_left)
                               : (trace_cell.current = trace_directions::diagonal, score_cell.current);

        if constexpr (align_local_t::value)
        {
            score_cell.current =
                (score_cell.current < 0) ? (trace_cell.current = trace_directions::none, 0) : score_cell.current;
        }
        // Check if this was the optimum. Possibly a noop.
        static_cast<alignment_algorithm_t const &>(*this).check_score_of_cell(current_cell, cache);

        // At the top of the band we can not come from up but only diagonal or left, so the next vertical must be a
        // gap open.
        score_cell.up = score_cell.current + cache.gap_open_score; // add gap open cost
        trace_cell.up = trace_directions::up_open;
    }

    /*!\brief Initialise the alignment state for affine gap computation.
     * \tparam alignment_configuration_t The type of alignment configuration.
     * \param[in] config The alignment configuration.
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
        auto scheme =
            config.get_or(align_cfg::gap_cost_affine{align_cfg::open_score{0}, align_cfg::extension_score{-1}});

        alignment_state.gap_extension_score = static_cast<score_t>(scheme.extension_score);
        alignment_state.gap_open_score = static_cast<score_t>(scheme.extension_score + scheme.open_score);
    }

    alignment_state_t alignment_state{}; //!< The internal alignment state tracking the current alignment optimum.
};

} // namespace seqan3::detail
