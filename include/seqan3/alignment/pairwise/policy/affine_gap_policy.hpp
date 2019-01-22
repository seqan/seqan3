// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::affine_gap_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// affine_gap_policy
// ----------------------------------------------------------------------------

/*!\brief The hot kernel implementation using affine gaps.
 * \ingroup alignment_policy
 * \tparam derived_t   The derived alignment algorithm.
 * \tparam cell_t      The cell type of the dynamic programming matrix.
 */
template <typename derived_t, typename cell_t>
class affine_gap_policy
{
protected:
    /*!\name Member types
     * \{
     */
    //!\brief The underlying score type.
    using score_t = std::tuple_element_t<0, cell_t>;
    //!\}

    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr affine_gap_policy()                                      noexcept = default;
    constexpr affine_gap_policy(affine_gap_policy const &)             noexcept = default;
    constexpr affine_gap_policy(affine_gap_policy &&)                  noexcept = default;
    constexpr affine_gap_policy & operator=(affine_gap_policy const &) noexcept = default;
    constexpr affine_gap_policy & operator=(affine_gap_policy &&)      noexcept = default;
    ~affine_gap_policy()                                               noexcept = default;
    //!}

    /*!\brief Computes the score for current cell.
     * \tparam        score_cell_type The type of the current cell.
     * \tparam        cache_t         The type of the cache.
     * \param[in,out] current_cell    The current cell in the dynamic programming matrix.
     * \param[in,out] cache           The cache storing hot helper variables.
     * \param[in]     score           The score of comparing the respective letters of the first and the second sequence.
     */
    template <typename score_cell_type, typename cache_t>
    constexpr void compute_cell(score_cell_type & current_cell,
                                cache_t & cache,
                                score_t const score) const noexcept
    {
        // Unpack the cached variables.
        auto & [main_score, hz_score] = current_cell;
        auto & [prev_cell, gap_open, gap_extend] = cache;
        auto & [prev_score, vt_score] = prev_cell;

        //TODO For the local alignment we might be able to use GCC overlow builtin arithmetics, which
        // allows us to check if overflow/underflow would happen. Not sure, if this helps with the performance though.
        // (see https://gcc.gnu.org/onlinedocs/gcc/Integer-Overflow-Builtins.html)
        score_t tmp = prev_score + score;
        tmp = std::max(tmp, vt_score);
        tmp = std::max(tmp, hz_score);
        prev_score = main_score;
        main_score = tmp;
        tmp += gap_open;
        vt_score += gap_extend;
        hz_score += gap_extend;
        vt_score = std::max(vt_score, tmp);
        hz_score = std::max(hz_score, tmp);
    }

    /*!\brief Creates the cache used for affine gap computation.
     * \tparam    gap_scheme_t The type of the gap scheme.
     * \param[in] scheme       The configured gap scheme.
     * \returns The created cache.
     */
    template <typename gap_scheme_t>
    constexpr auto make_cache(gap_scheme_t && scheme) const noexcept
    {
        return std::tuple{cell_t{},
                          static_cast<score_t>(scheme.get_gap_open_score() + scheme.get_gap_score()),
                          static_cast<score_t>(scheme.get_gap_score())};
    }
};

} // namespace seqan3::detail
