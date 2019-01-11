// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::unbanded_score_trace_dp_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alignment/pairwise/policy/unbanded_dp_matrix_policy.hpp>
#include <seqan3/std/span.hpp>

namespace seqan3::detail
{
/*!\brief Manages the allocation and provision of an unbanded dynamic programming matrix.
 * \ingroup alignment_policy
 * \tparam derived_t   The derived alignment algorithm.
 * \tparam score_allocator_t The allocator type used for allocating the score matrix.
 * \tparam trace_allocator_t The allocator type used for allocating the trace matrix.
 */
template <typename derived_t, typename score_allocator_t, typename trace_allocator_t>
class unbanded_score_trace_dp_matrix_policy :
    public unbanded_dp_matrix_policy<unbanded_score_trace_dp_matrix_policy<derived_t,
                                                                           score_allocator_t,
                                                                           trace_allocator_t>,
                                     score_allocator_t>
{
private:

    //!\brief The base type
    using base_t = unbanded_dp_matrix_policy<unbanded_score_trace_dp_matrix_policy<derived_t,
                                                                                   score_allocator_t,
                                                                                   trace_allocator_t>,
                                             score_allocator_t>;

    //!\brief Befriends the derived class to grant it access to the private members.
    friend derived_t;
    //!\brief Make dimension_first_batch visible in this class.
    using base_t::dimension_first_batch;
    //!\brief Make dimension_second_batch visible in this class.
    using base_t::dimension_second_batch;
    //!\brief Make score_matrix visible in this class.
    using base_t::score_matrix;

    /*!\name Member types
     * \{
     */
    //!\brief The underlying cell type of the dynamic programming matrix.
    using cell_type  = typename score_allocator_t::value_type;
    using trace_type = typename trace_allocator_t::value_type;
    //!\}

    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr unbanded_score_trace_dp_matrix_policy() = default;
    constexpr unbanded_score_trace_dp_matrix_policy(unbanded_score_trace_dp_matrix_policy const &) = default;
    constexpr unbanded_score_trace_dp_matrix_policy(unbanded_score_trace_dp_matrix_policy &&) = default;
    constexpr unbanded_score_trace_dp_matrix_policy & operator=(unbanded_score_trace_dp_matrix_policy const &) = default;
    constexpr unbanded_score_trace_dp_matrix_policy & operator=(unbanded_score_trace_dp_matrix_policy &&) = default;
    ~unbanded_score_trace_dp_matrix_policy() = default;
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
        dimension_first_batch = seqan3::size(std::forward<first_batch_t>(first_batch)) + 1;
        dimension_second_batch = seqan3::size(std::forward<second_batch_t>(second_batch)) + 1;

        // We use only one column to compute the score.
        score_matrix.resize(dimension_second_batch);
        // We use the full matrix to store the trace direction.
        trace_matrix.resize(dimension_first_batch * dimension_second_batch);
        trace_matrix_iter = trace_matrix.begin();
    }

    //!\brief Returns the active column that is computed next.
    constexpr auto active_column()
    {
        return std::tuple{std::span{score_matrix}, std::span{&(*trace_matrix_iter), dimension_second_batch}};
    }

    //!\brief Moves internal matrix pointer to the next column.
    constexpr void next_column() noexcept
    {
        trace_matrix_iter += dimension_second_batch;
    }

    //!\brief The data container.
    std::vector<trace_type, trace_allocator_t> trace_matrix{};
    //!\brief The current iterator in the trace matrix.
    typename std::vector<trace_type, trace_allocator_t>::iterator trace_matrix_iter{};
};
} // namespace seqan3::detail
