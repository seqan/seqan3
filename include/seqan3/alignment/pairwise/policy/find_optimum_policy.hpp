// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::find_optimum_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <range/v3/algorithm/for_each.hpp>

#include <seqan3/alignment/matrix/alignment_optimum.hpp>
#include <seqan3/alignment/pairwise/detail/alignment_algorithm_cache.hpp>
#include <seqan3/core/type_traits/deferred_crtp_base.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief The default traits class for seqan3::detail::find_optimum_policy.
 * \ingroup alignment_policy
 *
 * \details
 *
 * Defines the behaviour of a global alignment in which only the last cell of the dynamic programming matrix is
 * checked for the optimum.
 */
struct default_find_optimum_trait
{
    //!\brief Disables optimum search in every cell of the dynamic programming matrix.
    using find_in_every_cell_type  = std::false_type;
    //!\brief Disables optimum search in the last row of the dynamic programming matrix.
    using find_in_last_row_type    = std::false_type;
    //!\brief Disables optimum search in the last column of the dynamic programming matrix.
    using find_in_last_column_type = std::false_type;
};

/*!\brief The CRTP-policy to determine the optimum of the dynamic programming matrix.
 * \ingroup alignment_policy
 * \tparam alignment_algorithm_t The derived type (seqan3::detail::alignment_algorithm) to be augmented with this
 *                               CRTP-policy.
 * \tparam traits_type A traits type that determines which cells should be considered for the optimum.
 *                     Defaults to seqan3::detail::default_find_optimum_trait.
 *
 * \details
 *
 * This class determines the matrix wide optimum. The search space can be further refined using the
 * `traits_type` which configures the search space of the alignment matrix.
 */
template <typename derived_t, typename traits_type = default_find_optimum_trait>
class find_optimum_policy
{
private:
    //!\brief Befriends the derived class to grant it access to the private members.
    friend derived_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr find_optimum_policy() = default;                                        //!< Defaulted
    constexpr find_optimum_policy(find_optimum_policy const &) = default;             //!< Defaulted
    constexpr find_optimum_policy(find_optimum_policy &&) = default;                  //!< Defaulted
    constexpr find_optimum_policy & operator=(find_optimum_policy const &) = default; //!< Defaulted
    constexpr find_optimum_policy & operator=(find_optimum_policy &&) = default;      //!< Defaulted
    ~find_optimum_policy() = default;                                                 //!< Defaulted
    //!\}

protected:
    /*!\brief Checks every cell of the dynamic programming matrix.
     * \tparam score_t The type of the score.
     * \param[in] current_cell The current cell.
     * \param[in,out] cache The cache with the current optimum to update.
     *
     * \details
     *
     * This function resolves to a "NO-OP" function if the trait for searching every cell is set to std::false_type.
     */
    template <typename cell_t, typename score_t>
    constexpr void check_score([[maybe_unused]] cell_t const & current_cell,
                               [[maybe_unused]] alignment_algorithm_cache<score_t> & cache) const noexcept
    {
        using std::get;

        if constexpr (traits_type::find_in_every_cell_type::value)
        {
            cache.optimum = (get<0>(current_cell).current > cache.optimum.score)
                            ? alignment_optimum{get<0>(current_cell).current, get<1>(current_cell).coordinate}
                            : cache.optimum;
        }

    }

    /*!\brief Checks a cell of the last row of the dynamic programming matrix.
     * \tparam score_t The type of the score.
     * \param[in] current_cell The current cell.
     * \param[in,out] cache The cache with the current optimum to update.
     *
     * \details
     *
     * This function resolves to a "NO-OP" function if the trait for searching the last row is set to std::false_type.
     * Due to a column based iteration layout this computes only one cell at a time. The alignment algorithm
     * takes care of calling this function for the appropriate cells.
     */
    template <typename cell_t, typename score_t>
    constexpr void check_score_last_row([[maybe_unused]] cell_t const & current_cell,
                                        [[maybe_unused]] alignment_algorithm_cache<score_t> & cache) const noexcept
    {
        using std::get;

        if constexpr (traits_type::find_in_last_row_type::value)
        {
            cache.optimum = (get<0>(current_cell).current > cache.optimum.score)
                            ? alignment_optimum{get<0>(current_cell).current, get<1>(current_cell).coordinate}
                            : cache.optimum;
        }
    }

    /*!\brief Checks the complete last column for the optimal score.
     * \tparam rng_t   The type of the last column; must model std::ranges::bidirectional_range.
     * \tparam score_t The type of the optimal score.
     * \param[in] current_cell The current cell.
     * \param[in,out] cache The cache with the current optimum to update.
     *
     * \details
     *
     * This function checks only the last element of the column (the score for the global alignment)
     * if the trait for searching the last column is set to std::false_type.
     * Due to a column based iteration layout the entire last column can be searched at once.
     */
    template <typename cell_t, typename score_t>
    constexpr void check_score_last_column_or_cell([[maybe_unused]] cell_t const & current_cell,
                                                   alignment_algorithm_cache<score_t> & cache) const noexcept
    {
        using std::get;
        // Only check the entire column if it was configured to search here.
        if constexpr (traits_type::find_in_last_column_type::value)
        {
            // get the iterator from the class itself.
            derived_t const * me = static_cast<derived_t const *>(this);
            std::ranges::for_each(views::zip(*me->score_matrix_iter, *me->trace_matrix_iter), [&](auto && cell)
            {
                cache.optimum = (get<0>(cell).current > cache.optimum.score)
                                ? alignment_optimum{get<0>(cell).current, get<1>(cell).coordinate}
                                : cache.optimum;
            });
        }
        else  // Only check the last cell for the global alignment.
        {
            cache.optimum = (get<0>(current_cell).current > cache.optimum.score)
                                ? alignment_optimum{get<0>(current_cell).current, get<1>(current_cell).coordinate}
                                : cache.optimum;
        }
    }
};

} // namespace seqan3::detail
