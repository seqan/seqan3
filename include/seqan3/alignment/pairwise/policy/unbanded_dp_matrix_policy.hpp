// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::unbanded_dp_matrix_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <vector>
#include <tuple>

#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/std/span.hpp>

namespace seqan3::detail
{

/*!\brief Manages the allocation and provision of an unbanded dynamic programming matrix.
 * \ingroup alignment_policy
 * \tparam derived_t   The derived alignment algorithm.
 * \tparam allocator_t The allocator type used for allocating the dynamic programming matrix.
 */
template <typename derived_t, typename allocator_t>
class unbanded_dp_matrix_policy
{
private:

    //!\brief Befriends the derived class to grant it access to the private members.
    friend derived_t;

    /*!\name Member types
     * \{
     */
    //!\brief The underlying cell type of the dynamic programming matrix.
    using cell_type = typename allocator_t::value_type;
    //!\}

    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr unbanded_dp_matrix_policy()                                              = default;
    constexpr unbanded_dp_matrix_policy(unbanded_dp_matrix_policy const &)             = default;
    constexpr unbanded_dp_matrix_policy(unbanded_dp_matrix_policy &&)                  = default;
    constexpr unbanded_dp_matrix_policy & operator=(unbanded_dp_matrix_policy const &) = default;
    constexpr unbanded_dp_matrix_policy & operator=(unbanded_dp_matrix_policy &&)      = default;
    ~unbanded_dp_matrix_policy()                                                       = default;
    //!}

    /*!\brief Allocates the memory for the dynamic programming matrix given the two sequences.
     * \tparam first_batch_t   The type of the first sequence (or packed sequences).
     * \tparam second_batch_t  The type of the second sequence (or packed sequences).
     * \param[in] first_batch  The first sequence (or packed sequences).
     * \param[in] second_batch The first sequence (or packed sequences).
     */
    template <typename first_batch_t, typename second_batch_t>
    constexpr void allocate_matrix(first_batch_t && first_batch, second_batch_t && second_batch)
    {
        dimension_first_batch = std::ranges::size(std::forward<first_batch_t>(first_batch)) + 1;
        dimension_second_batch = std::ranges::size(std::forward<second_batch_t>(second_batch)) + 1;

        // We use only one column to compute the score.
        score_matrix.resize(dimension_second_batch);
    }

    //!\brief Returns the active column that is computed next.
    constexpr auto active_column()
    {
        return std::tuple{std::span{score_matrix}, std::ignore};
    }

    //!\brief Moves internal matrix pointer to the next column.
    constexpr void next_column()
    {
        // noop
    }

    //!\brief The data container.
    std::vector<cell_type, allocator_t> score_matrix{};
    //!\brief Caches the size of the horizontal dimension (number of columns).
    size_t dimension_first_batch  = 0;
    //!\brief Caches the size of the vertical dimension (number of rows).
    size_t dimension_second_batch = 0;
};
} // namespace seqan3::detail
