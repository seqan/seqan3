// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
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

#include <seqan3/core/metafunction/deferred_crtp_base.hpp>
#include <seqan3/range/shortcuts.hpp>

namespace seqan3::detail
{

/*!\brief The default traits class for seqan3::detail::find_optimum_policy.
 * \ingroup alignment_policy
 *
 * \details
 *
 * Defines the behavior of a global alignment, in which only the last cell of the dynamic programming matrix is
 * checked for the optimum.
 */
struct default_find_optimum_trait
{
    using find_in_every_cell_t  = std::false_type;
    using find_in_last_row_t    = std::false_type;
    using find_in_last_column_t = std::false_type;
};

/*!\brief A policy to determine the optimum of the dynamic programming matrix.
 * \ingroup alignment_policy
 * \tparam derived_t        The type of the derived class.
 * \tparam core_policy_t    The type of the core alignment policy which also might need access to the member functions.
 * \tparam traits_type      A traits type that determines which cells should be considered for the optimum.
 *                          Defaults to seqan3::detail::default_find_optimum_trait.
 */
template <typename derived_t, typename core_policy_t, typename traits_type = default_find_optimum_trait>
class find_optimum_policy
{

private:

    //!\brief Befriends the derived class to grant it access to the private members.
    friend derived_t;
    //!\brief Befriends the core_policy_t which is still a seqan3::detail::deferred_crtp_base and must be resolved
    //!\      accordingly.
    friend invoke_deferred_crtp_base<core_policy_t, derived_t>;

    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr find_optimum_policy()                                        = default;
    constexpr find_optimum_policy(find_optimum_policy const &)             = default;
    constexpr find_optimum_policy(find_optimum_policy &&)                  = default;
    constexpr find_optimum_policy & operator=(find_optimum_policy const &) = default;
    constexpr find_optimum_policy & operator=(find_optimum_policy &&)      = default;
    ~find_optimum_policy()                                                 = default;
    //!}

    /*!\brief Checks every cell of the dynamic programming matrix.
     * \tparam score_t  The type of the score.
     * \param[in]     val       The current value.
     * \param[in,out] optimum   The current optimum to compare with.
     */
    template <typename score_t>
    constexpr void check_score([[maybe_unused]] score_t const val,
                               [[maybe_unused]] score_t & optimum) const noexcept
    {
        if constexpr (traits_type::find_in_every_cell_t::value)
            optimum = std::max(optimum, val);
    }

    /*!\brief Checks the last cell of the last row of the dynamic programming matrix.
     * \tparam score_t  The type of the score.
     * \param[in]     val       The current value.
     * \param[in,out] optimum   The current optimum to compare with.
     */
    template <typename score_t>
    constexpr void check_score_last_row([[maybe_unused]] score_t const val,
                                        [[maybe_unused]] score_t & optimum) const noexcept
    {
        if constexpr (traits_type::find_in_last_row_t::value)
            optimum = std::max(optimum, val);
    }

    /*!\brief Checks the complete last column for the optimal score.
     * \tparam rng_t    The type of the last column. Must be at least an input range.
     * \tparam score_t  The type of the optimal score.
     * \param[in]     rng       The last column of the dynamic programming matrix.
     * \param[in,out] optimum   The current optimum to compare with.
     */
    template <typename rng_t, typename score_t>
    constexpr void check_score_last_column(rng_t const & rng, score_t & optimum) const noexcept
    {
        using std::get;
        // Only check the entire column if it was configured to search here.
        if constexpr (traits_type::find_in_last_column_t::value)
        {
            ranges::for_each(rng, [&](auto && tpl)
            {
                optimum = std::max(optimum, get<0>(tpl));
            });
        }
        else  // Only check the last cell of the global alignment.
        {
            auto val = get<0>(*(seqan3::end(rng) - 1));
            optimum = std::max(optimum, val);
        }
    }
};

} // namespace seqan3::detail
