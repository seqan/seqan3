// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::find_optimum_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <range/v3/algorithm/for_each.hpp>

#include <seqan3/alignment/matrix/alignment_optimum.hpp>
#include <seqan3/core/metafunction/deferred_crtp_base.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/std/concepts>
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
    //!\brief Disables lookup of the optimal score in every cell.
    using find_in_every_cell_type  = std::false_type;
    //!\brief Disables lookup of the optimal score in the last row.
    using find_in_last_row_type    = std::false_type;
    //!\brief Disables lookup of the optimal score in the last column.
    using find_in_last_column_type = std::false_type;
};

/*!\brief A policy to determine the optimum of the dynamic programming matrix.
 * \ingroup alignment_policy
 * \tparam derived_t        The type of the derived class.
 * \tparam traits_type      A traits type that determines which cells should be considered for the optimum.
 *                          Defaults to seqan3::detail::default_find_optimum_trait.
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

    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr find_optimum_policy() = default;
    constexpr find_optimum_policy(find_optimum_policy const &) = default;
    constexpr find_optimum_policy(find_optimum_policy &&) = default;
    constexpr find_optimum_policy & operator=(find_optimum_policy const &) = default;
    constexpr find_optimum_policy & operator=(find_optimum_policy &&) = default;
    ~find_optimum_policy() = default;
    //!\}

protected:

    /*!\brief Checks every cell of the dynamic programming matrix.
     * \tparam score_t The type of the score.
     * \param[in]     current   The current value.
     * \param[in,out] optimum   The current optimum to compare with.
     *
     * \details
     *
     * This function resolves to a "NO-OP" function if the trait for searching every cell is set to std::false_type.
     */
    template <typename score_t>
    constexpr void check_score([[maybe_unused]] alignment_optimum<score_t> const & current,
                               [[maybe_unused]] alignment_optimum<score_t> & optimum) const noexcept
    {
        if constexpr (traits_type::find_in_every_cell_type::value)
        {
            optimum = std::max(current, optimum, alignment_optimum_compare_less{});
        }

    }

    /*!\brief Checks a cell of the last row of the dynamic programming matrix.
     * \tparam score_t The type of the score.
     * \param[in]     current   The current value.
     * \param[in,out] optimum   The current optimum to compare with.
     *
     * \details
     *
     * This function resolves to a "NO-OP" function if the trait for searching the last row is set to std::false_type.
     * Due to a column based iteration layout this computes only one cell at a time. The alignment algorithm
     * takes care of calling this function for the appropriate cells.
     */
    template <typename score_t>
    constexpr void check_score_last_row([[maybe_unused]] alignment_optimum<score_t> const & current,
                                        [[maybe_unused]] alignment_optimum<score_t> & optimum) const noexcept
    {
        if constexpr (traits_type::find_in_last_row_type::value)
            optimum = std::max(current, optimum, alignment_optimum_compare_less{});
    }

    /*!\brief Checks the complete last column for the optimal score.
     * \tparam rng_t   The type of the last column; must model std::ranges::BidirectionalRange.
     * \tparam score_t The type of the optimal score.
     * \param[in]     rng       The last column of the dynamic programming matrix.
     * \param[in,out] optimum   The current optimum to compare with.
     *
     * \details
     *
     * This function checks only the last element of the column (the score for the global alignment)
     * if the trait for searching the last column is set to std::false_type.
     * Due to a column based iteration layout the entire last column can be searched at once.
     */
    template <std::ranges::BidirectionalRange rng_t, typename score_t>
    constexpr void check_score_last_column(rng_t const & rng,
                                           alignment_optimum<score_t> & optimum) const noexcept
    {
        using std::get;
        // Only check the entire column if it was configured to search here.
        if constexpr (traits_type::find_in_last_column_type::value)
        {
            ranges::for_each(rng, [&](auto && tpl)
            {
                optimum = std::max(alignment_optimum<score_t>{get<0>(get<0>(tpl)),
                                                              static_cast<alignment_coordinate>(get<1>(tpl))},
                                   optimum,
                                   alignment_optimum_compare_less{});
            });
        }
        else  // Only check the last cell for the global alignment.
        {
            auto && current = *(--seqan3::end(rng));
            optimum = std::max(alignment_optimum<score_t>{get<0>(get<0>(current)),
                                                          static_cast<alignment_coordinate>(get<1>(current))},
                               optimum,
                               alignment_optimum_compare_less{});
        }
    }
};

} // namespace seqan3::detail
