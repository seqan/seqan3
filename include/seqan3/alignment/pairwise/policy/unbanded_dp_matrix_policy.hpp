// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::unbanded_dp_matrix_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <vector>
#include <utility>

#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/std/view/subrange.hpp>
#include <seqan3/std/ranges>

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
protected:

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
    void allocate_score_matrix(first_batch_t && first_batch, second_batch_t && second_batch)
    {
        dimension_first_batch = std::ranges::size(first_batch);
        dimension_second_batch = std::ranges::size(second_batch);

        // We use only one column to compute the score.
        score_matrix.resize(dimension_second_batch + 1);
    }

    //!\brief Returns the current column of the alignment matrix.
    auto current_column()
    {
        using iter_t = std::ranges::iterator_t<decltype(score_matrix)>;
        return view::subrange<iter_t, iter_t>{std::ranges::begin(score_matrix), std::ranges::end(score_matrix)};
    }

    //!\brief The data container.
    std::vector<cell_type, allocator_t> score_matrix{};
    //!\brief Caches the size of the horizontal dimension (number of columns).
    size_t dimension_first_batch  = 0;
    //!\brief Caches the size of the vertical dimension (number of rows).
    size_t dimension_second_batch = 0;
};
} // namespace seqan3::detail
