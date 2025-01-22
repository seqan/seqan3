// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::policy_affine_gap_recursion.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/alignment/matrix/detail/affine_cell_proxy.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/concept.hpp>

namespace seqan3::detail
{

/*!\brief Implements the alignment recursion function for the alignment algorithm using affine gap costs.
 * \ingroup alignment_pairwise
 *
 * \tparam alignment_configuration_t The type of the alignment configuration.
 *
 * \details
 *
 * Implements the functions to initialise and compute the alignment matrix using the recursion formula for affine gaps.
 * Other policies can inherit from this policy and overload the recursion functions, e.g. to change the
 * initialisation of the alignment matrix.
 *
 * \note For more information, please refer to the original article for the alignment with affine gap cost function:
 *       GOTOH, Osamu. An improved algorithm for matching biological sequences.
 *       Journal of molecular biology, 1982, 162. Jg., Nr. 3, S. 705-708.
 */
template <typename alignment_configuration_t>
class policy_affine_gap_recursion
{
protected:
    //!\brief The configuration traits type.
    using traits_type = alignment_configuration_traits<alignment_configuration_t>;
    //!\brief The configured original score type.
    using original_score_type = typename traits_type::original_score_type;
    //!\brief The configured score type.
    using score_type = typename traits_type::score_type;
    //!\brief The internal tuple storing the scores of an affine cell.
    using affine_score_tuple_t = std::tuple<score_type, score_type, score_type>;
    //!\brief The affine cell type returned by the functions.
    using affine_cell_type = affine_cell_proxy<affine_score_tuple_t>;

    //!\brief The score for a gap extension.
    score_type gap_extension_score{};
    //!\brief The score for a gap opening including the gap extension.
    score_type gap_open_score{};

    //!\brief Initialisation state of the first row of the alignment.
    bool first_row_is_free{};
    //!\brief Initialisation state of the first column of the alignment.
    bool first_column_is_free{};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_affine_gap_recursion() = default;                                                //!< Defaulted.
    policy_affine_gap_recursion(policy_affine_gap_recursion const &) = default;             //!< Defaulted.
    policy_affine_gap_recursion(policy_affine_gap_recursion &&) = default;                  //!< Defaulted.
    policy_affine_gap_recursion & operator=(policy_affine_gap_recursion const &) = default; //!< Defaulted.
    policy_affine_gap_recursion & operator=(policy_affine_gap_recursion &&) = default;      //!< Defaulted.
    ~policy_affine_gap_recursion() = default;                                               //!< Defaulted.

    /*!\brief Construction and initialisation using the alignment configuration.
     * \param[in] config The alignment configuration.
     *
     * \details
     *
     * Initialises the gap open score and gap extension score for this policy.
     * If no gap cost model was provided by the user the default gap costs `-10` and `-1` are set for the gap open score
     * and the gap extension score respectively.
     */
    explicit policy_affine_gap_recursion(alignment_configuration_t const & config)
    {
        // Get the gap scheme from the config or choose -1 and -10 as default.
        auto const & selected_gap_scheme =
            config.get_or(align_cfg::gap_cost_affine{align_cfg::open_score{-10}, align_cfg::extension_score{-1}});

        gap_extension_score = maybe_convert_to_simd(selected_gap_scheme.extension_score);
        gap_open_score = maybe_convert_to_simd(selected_gap_scheme.open_score) + gap_extension_score;

        auto method_global_config = config.get_or(align_cfg::method_global{});
        first_row_is_free = method_global_config.free_end_gaps_sequence1_leading;
        first_column_is_free = method_global_config.free_end_gaps_sequence2_leading;
    }
    //!\}

    /*!\brief Computes an inner cell of the alignment matrix.
     *
     * \tparam affine_cell_t The type of the affine cell; must be an instance of seqan3::detail::affine_cell_proxy.
     *
     * \param[in] diagonal_score The previous diagonal score, which corresponds to \f$M[i - 1, j - 1]\f$.
     * \param[in] previous_cell The predecessor cell corresponding to the values \f$V[i - 1, j]\f$ and \f$H[i, j -1]\f$.
     * \param[in] sequence_score The score obtained from the scoring scheme for the current cell (\f$ \delta\f$).
     *
     * \returns The computed affine cell.
     *
     * \details
     *
     * Computes the current cell according to following recursion formula:
     * * \f$ H[i, j] = \max \{M[i, j - 1] + g_o, H[i, j - 1] + g_e\}\f$
     * * \f$ V[i, j] = \max \{M[i - 1, j] + g_o, V[i - 1, j] + g_e\}\f$
     * * \f$ M[i, j] = \max \{M[i - 1, j - 1] + \delta, H[i, j], V[i, j]\}\f$
     */
    template <typename affine_cell_t>
    affine_cell_type compute_inner_cell(score_type diagonal_score,
                                        affine_cell_t previous_cell,
                                        score_type const sequence_score) const noexcept
    {
        diagonal_score += sequence_score;
        score_type horizontal_score = previous_cell.horizontal_score();
        score_type vertical_score = previous_cell.vertical_score();

        diagonal_score = (diagonal_score < vertical_score) ? vertical_score : diagonal_score;
        diagonal_score = (diagonal_score < horizontal_score) ? horizontal_score : diagonal_score;

        score_type tmp = diagonal_score + gap_open_score;
        vertical_score += gap_extension_score;
        horizontal_score += gap_extension_score;

        // store the vertical_score and horizontal_score value in the next path
        vertical_score = (vertical_score < tmp) ? tmp : vertical_score;
        horizontal_score = (horizontal_score < tmp) ? tmp : horizontal_score;

        return {diagonal_score, horizontal_score, vertical_score};
    }

    /*!\brief Initialises the first cell of the alignment matrix in the top left corner of the matrix.
     *
     * \returns The computed affine cell.
     *
     * \details
     *
     * Initialises the cell at the origin of the alignment matrix (top left corner of the matrix). The optimal score is
     * initialised to 0, while the value of the horizontal and vertical matrix are initialised as:
     * \f$V[0, 0] = H[0, 0] = g_o\f$.
     */
    affine_cell_type initialise_origin_cell() const noexcept
    {
        return {score_type{},
                first_row_is_free ? score_type{} : gap_open_score,
                first_column_is_free ? score_type{} : gap_open_score};
    }

    /*!\brief Initialises a cell of the first alignment matrix column.
     *
     * \tparam affine_cell_t The type of the affine cell; must be an instance of seqan3::detail::affine_cell_proxy.
     *
     * \param[in] previous_cell The predecessor cell on the same column \f$M[i-1, 0]\f$.
     *
     * \returns The computed affine cell.
     *
     * \details
     *
     * Initialises a cell of the first alignment matrix column. The optimal score is the same as the vertical score
     * which is equal to \f$V[i, 0] = M[i, 0] = g_o + g_e * i\f$. The horizontal score is initialised to
     * \f$H[i, 0] = V[i, 0] + g_o\f$ to prohibit extending a gap in the horizontal matrix from \f$H[i, 0]\f$.
     */
    template <typename affine_cell_t>
    affine_cell_type initialise_first_column_cell(affine_cell_t previous_cell) const noexcept
    {
        score_type new_vertical = previous_cell.vertical_score() + gap_extension_score;
        return {previous_cell.vertical_score(),
                previous_cell.vertical_score() + gap_open_score,
                first_column_is_free ? previous_cell.vertical_score() : new_vertical};
    }

    /*!\brief Initialises the first cell of a alignment matrix column.
     *
     * \tparam affine_cell_t The type of the affine cell; must be an instance of seqan3::detail::affine_cell_proxy.
     *
     * \param[in] previous_cell The predecessor cell on the same row \f$M[0, j-1]\f$.
     *
     * \returns The computed affine cell.
     *
     * \details
     *
     * Initialises the first cell of a alignment matrix column. The optimal score is the same as the horizontal score
     * which is equal to \f$H[0, j] = M[0, j] = g_o + g_e * j\f$. The vertical score is initialised to
     * \f$V[0,j] = H[0, j] + g_o\f$ to prohibit extending a gap in the vertical matrix from \f$V[0, j]\f$.
     */
    template <typename affine_cell_t>
    affine_cell_type initialise_first_row_cell(affine_cell_t previous_cell) const noexcept
    {
        score_type new_horizontal_score = previous_cell.horizontal_score() + gap_extension_score;
        return {previous_cell.horizontal_score(),
                first_row_is_free ? previous_cell.horizontal_score() : new_horizontal_score,
                previous_cell.horizontal_score() + gap_open_score};
    }

    /*!\brief Returns the lowest viable score.
     *
     * \details
     *
     * In some versions of the algorithms a value representing minus infinity is needed. Since the data type is an
     * signed integral there is no infinity but only the lowest possible value that can be represented by the score
     * type. In order to avoid unnecessary if conditions to protect against signed integer underflow the lowest viable
     * score is computed. Subtracting a gap penalty from this will still result in a valid score which represents
     * minus infinity.
     */
    score_type lowest_viable_score() const noexcept
    {
        if constexpr (simd_concept<score_type>)
            assert(gap_open_score[0] <= 0 && gap_extension_score[0] <= 0);
        else
            assert(gap_open_score <= 0 && gap_extension_score <= 0);

        return maybe_convert_to_simd(std::numeric_limits<original_score_type>::lowest())
             - (gap_open_score + gap_extension_score);
    }

    /*!\brief Converts the given score type to a simd vector if the alignment is executed in vectorised mode.
     *
     * \tparam score_t The score type to convert; must model seqan3::arithmetic.
     * \param[in] score The score to convert.
     *
     * \returns The score converted to the target simd vector or the unmodified value if in scalar mode.
     */
    template <typename score_t>
        requires arithmetic<std::remove_cvref_t<score_t>>
    constexpr auto maybe_convert_to_simd(score_t && score) const noexcept
    {
        if constexpr (simd_concept<score_type>)
            return simd::fill<score_type>(std::forward<score_t>(score));
        else // Return unmodified.
            return std::forward<score_t>(score);
    }
};
} // namespace seqan3::detail
