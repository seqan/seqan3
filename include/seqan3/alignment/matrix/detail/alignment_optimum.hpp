// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::alignment_optimum.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/concept.hpp>
#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>

namespace seqan3::detail
{

/*!\brief Stores the current optimum of the alignment algorithm.
 * \ingroup alignment_matrix
 *
 * \tparam score_t The type of the tracked alignment score; must model either seqan3::arithmetic or
 *                 seqan3::simd_concept.
 *
 * \details
 *
 * Stores the optimal score of the alignment computation and the corresponding indices of the cell with the optimal
 * score within the alignment matrix.
 * In case the optimum is used for the vectorised alignment computation this optimum stores the optimal scores and the
 * respective cells as simd vectors.
 */
template <typename score_t>
struct alignment_optimum
#if SEQAN3_DOXYGEN_ONLY(1) 0
{
    //!\brief The index type used to store the alignment coordinates of the optimum.
    using index_t = IMPLEMENTATION_DEFINED;

    //!\brief The index of the alignment matrix column.
    index_t column_index{};
    //!\brief The index of the alignment matrix row.
    index_t row_index{};
    //!\brief The optimal score whose initialisation is implementation defined.
    score_t score = IMPLEMENTATION_DEFINED;

    /*!\brief Compares the score with the given score and updates the optimum if the new score is bigger than
     *        the current one.
     *
     * \tparam column_index_t The index type for the column index; must model std::unsigned_integral.
     * \tparam row_index_t The index type for the row index; must model std::unsigned_integral.
     *
     * \param[in] compare_score The new score to compare with.
     * \param[in] column_index The respective column index of the alignment matrix.
     * \param[in] row_index The respective row index of the alignment matrix.
     *
     * \details
     *
     * Only updates the current optimum if the new score is greater than the current one. Note in the case of computing
     * a vectorised alignment only the positions of the simd vector are updated whose score is greater than the current
     * scores.
     */
    template <typename column_index_t, typename row_index_t>
    void update_if_new_optimal_score(score_t const & compare_score,
                                     column_index_type<column_index_t> column_index,
                                     row_index_type<row_index_t> row_index) noexcept;
}
#endif //SEQAN3_DOXYGEN_ONLY(1): This code block is only dis
;

//!\cond
template <arithmetic score_t>
struct alignment_optimum<score_t>
{
    size_t column_index{};
    size_t row_index{};
    score_t score{std::numeric_limits<score_t>::lowest()};

    template <std::integral column_index_t, std::integral row_index_t>
    constexpr void update_if_new_optimal_score(score_t const & compare_score,
                                               column_index_type<column_index_t> column_index,
                                               row_index_type<row_index_t> row_index) noexcept
    {
        score = (compare_score > score)
                  ? (this->column_index = column_index.get(), this->row_index = row_index.get(), compare_score)
                  : score;
    }
};

template <simd_concept score_t>
struct alignment_optimum<score_t>
{
    using scalar_t = typename simd_traits<score_t>::scalar_type;

    score_t column_index{};
    score_t row_index{};
    score_t score{simd::fill<score_t>(std::numeric_limits<scalar_t>::lowest())};

    template <std::integral column_index_t, std::integral row_index_t>
    constexpr void update_if_new_optimal_score(score_t const & compare_score,
                                               column_index_type<column_index_t> column_index,
                                               row_index_type<row_index_t> row_index) noexcept
    {
        auto mask = compare_score > score;
        score = mask ? compare_score : score;
        this->column_index = mask ? simd::fill<score_t>(column_index.get()) : this->column_index;
        this->row_index = mask ? simd::fill<score_t>(row_index.get()) : this->row_index;
    }
};
//!\endcond

/*!\name Type deduction guides
 * \relates seqan3::detail::alignment_optimum
 * \{
 */
//!\brief Default constructed objects deduce to `int32_t`.
alignment_optimum() -> alignment_optimum<int32_t>;

//!\brief Construction from column index, row index and the score deduces the score type.
template <typename column_index_t, typename row_index_t, typename score_t>
alignment_optimum(column_index_t, row_index_t, score_t) -> alignment_optimum<score_t>;
//!\}

} // namespace seqan3::detail
