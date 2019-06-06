// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_optimum.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

/*!\brief Stores the current optimum of the alignment algorithms.
 * \ingroup alignment_matrix
 * \tparam score_t The type of the tracked alignment score.
 *
 * \details
 *
 * This is an aggregate type, so the score needs to be passed before the seqan3::alignment_coordinate during
 * construction.
 */
template <Arithmetic score_t>
struct alignment_optimum
{
    //!\brief The optimal score.
    score_t score{std::numeric_limits<score_t>::lowest()};
    //!\brief The corresponding coordinate within the alignment matrix.
    alignment_coordinate coordinate{};
};

/*!\name Type deduction guides
 * \{
 */
//!\brief Default constructed objects deduce to `int32_t`.
alignment_optimum() -> alignment_optimum<int32_t>;

//!\brief Deduce the score type.
template <Arithmetic score_t>
alignment_optimum(score_t const, alignment_coordinate const) ->
    alignment_optimum<std::remove_reference_t<score_t>>;
//!\}

/*!\brief A less than comparator for two seqan3::detail::alignment_optimum objects.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * This function object is used in std::max functions to compare two seqan3::detail::alignment_optimum objects.
 */
struct alignment_optimum_compare_less
{

    /*!\brief Function call operator that implements less than comparison.
     * \tparam lhs_t The type of the left-hand side operand. Must be a type specialisation of
     *               seqan3::detail::alignment_optimum.
     * \tparam rhs_t The type of the right-hand side operand. Must be a type specialisation of
     *               seqan3::detail::alignment_optimum.
     *
     * \param[in] lhs The left-hand side operand.
     * \param[in] rhs The right-hand side operand.
     *
     * \returns bool `true` if `lhs.score < rhs.score`, otherwise `false`.
     */
    template <typename lhs_t, typename rhs_t>
    //!\cond
        requires (is_type_specialisation_of_v<lhs_t, alignment_optimum> &&
                  is_type_specialisation_of_v<rhs_t, alignment_optimum>)
    //!\endcond
    constexpr bool operator()(lhs_t const & lhs, rhs_t const & rhs) const
    {
        return lhs.score < rhs.score;
    }
};

} // namespace seqan3::detail
