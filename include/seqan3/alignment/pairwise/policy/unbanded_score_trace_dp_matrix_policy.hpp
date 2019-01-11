// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::unbanded_score_trace_dp_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <range/v3/view/zip.hpp>

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
    //!\brief Make dimension_first_range visible in this class.
    using base_t::dimension_first_range;
    //!\brief Make dimension_second_range visible in this class.
    using base_t::dimension_second_range;
    //!\brief Make score_matrix visible in this class.
    using base_t::score_matrix;

    /*!\name Member types
     * \{
     */
    //!\brief The underlying cell type of the dynamic programming matrix.
    //!\brief The underlying cell type of the scoring matrix.
    using cell_type = typename score_allocator_t::value_type;
    //!\brief The underlying cell type of the trace matrix.
    using trace_type = typename trace_allocator_t::value_type;
    //!\brief The type of the score matrix.
    using score_matrix_type = std::vector<cell_type, score_allocator_t>;
    //!\brief The type of the trace matrix.
    using trace_matrix_type = std::vector<trace_type, trace_allocator_t>;
    //!\}

    /*!\name Constructor, destructor and assignment
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
     * \tparam first_range_t   The type of the first sequence (or packed sequences).
     * \tparam second_range_t  The type of the second sequence (or packed sequences).
     * \param[in] first_range  The first sequence (or packed sequences).
     * \param[in] second_range The second sequence (or packed sequences).
     */
    template <typename first_range_t, typename second_range_t>
    constexpr void allocate_matrix(first_range_t && first_range, second_range_t && second_range)
    {
        dimension_first_range = seqan3::size(std::forward<first_range_t>(first_range)) + 1;
        dimension_second_range = seqan3::size(std::forward<second_range_t>(second_range)) + 1;

        // We use only one column to compute the score.
        score_matrix.resize(dimension_second_range);
        // We use the full matrix to store the trace direction.
        trace_matrix.resize(dimension_first_range * dimension_second_range);
        trace_matrix_iter = trace_matrix.begin();
    }

    //!\brief Returns the current column of the alignment matrix.
    constexpr auto current_column()
    {
        return ranges::view::zip(std::span{score_matrix}, std::span{&(*trace_matrix_iter), dimension_second_range});
    }

    //!\brief Moves internal matrix pointer to the next column.
    constexpr void next_column() noexcept
    {
        trace_matrix_iter += dimension_second_range;
    }

    //!\brief The data container.
    trace_matrix_type trace_matrix{};
    //!\brief The current iterator in the trace matrix.
    typename trace_matrix_type::iterator trace_matrix_iter{};
};
} // namespace seqan3::detail
