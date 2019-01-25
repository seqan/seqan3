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

#include <range/v3/view/iota.hpp>
#include <range/v3/view/zip.hpp>

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/pairwise/policy/banded_score_dp_matrix_policy.hpp>
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
class banded_score_trace_dp_matrix_policy :
    public banded_score_dp_matrix_policy<banded_score_trace_dp_matrix_policy<derived_t,
                                                                             score_allocator_t,
                                                                             trace_allocator_t>,
                                         score_allocator_t>
{
private:

    //!\brief The base type
    using base_t = banded_score_dp_matrix_policy<banded_score_trace_dp_matrix_policy<derived_t,
                                                                                     score_allocator_t,
                                                                                     trace_allocator_t>,
                                                 score_allocator_t>;

    //!\brief Befriends the derived class to grant it access to the private members.
    friend derived_t;

    // Import members from base class.
    using base_t::score_matrix;
    using base_t::dimension_first_range;
    using base_t::dimension_second_range;
    using base_t::current_column_index;
    using base_t::current_matrix_iter;
    using base_t::band_column_index;
    using base_t::band_row_index;

    // Import member functions from base class
    using base_t::current_band_size;
    using base_t::second_range_begin_offset;
    using base_t::band_touches_last_row;
    using base_t::trim_sequences;

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
    constexpr banded_score_trace_dp_matrix_policy() = default;
    constexpr banded_score_trace_dp_matrix_policy(banded_score_trace_dp_matrix_policy const &) = default;
    constexpr banded_score_trace_dp_matrix_policy(banded_score_trace_dp_matrix_policy &&) = default;
    constexpr banded_score_trace_dp_matrix_policy & operator=(banded_score_trace_dp_matrix_policy const &) = default;
    constexpr banded_score_trace_dp_matrix_policy & operator=(banded_score_trace_dp_matrix_policy &&) = default;
    ~banded_score_trace_dp_matrix_policy() = default;
    //!}

    /*!\brief Allocates the memory for the dynamic programming matrix given the two sequences.
     * \tparam first_range_t   The type of the first sequence (or packed sequences).
     * \tparam second_range_t  The type of the second sequence (or packed sequences).
     * \param[in] first_range  The first sequence (or packed sequences).
     * \param[in] second_range The second sequence (or packed sequences).
     */
    template <typename first_range_t, typename second_range_t, typename band_t>
    constexpr void allocate_matrix(first_range_t && first_range, second_range_t && second_range, band_t const & band)
    {
        base_t::allocate_matrix(first_range, second_range, band);

        // Allocate the k-banded score matrix for the trace matrix and set the matrix iterator to the first position
        // within the matrix.
        band_dimension = band_column_index + band_row_index + 1;
        trace_matrix.resize(band_dimension * dimension_first_range);
        trace_matrix_iter = seqan3::begin(trace_matrix) + band_column_index;
    }

    //!\brief Returns the current column of the alignment matrix.
    constexpr auto current_column() noexcept
    {
        auto span = base_t::current_band_size();

        assert(span > 0u);  // The span must always be greater than 0.

        // The begin coordinate in the current column begins at it - begin(matrix).
        // The end coordinate ends at it - begin(matrix) + current_band_size
        advanceable_alignment_coordinate<size_t, advanceable_alignment_coordinate_state::row>
            col_begin{column_index_type{current_column_index},
                      row_index_type{static_cast<size_t>(std::ranges::distance(seqan3::begin(score_matrix),
                                                                               current_matrix_iter))}};
        advanceable_alignment_coordinate<size_t, advanceable_alignment_coordinate_state::row>
            col_end{column_index_type{current_column_index}, row_index_type{col_begin.second_seq_pos + span}};

        // Return zip view over current column and current column shifted by one to access the previous horizontal.
        auto zip_score = ranges::view::zip(std::span{std::addressof(*current_matrix_iter), span},
                                           std::span{std::addressof(*(current_matrix_iter + 1)), span});
        return ranges::view::zip(std::move(zip_score),
                                 ranges::view::iota(col_begin, col_end),
                                 std::span{std::addressof(*trace_matrix_iter), span});
    }

    //!\brief Moves internal matrix pointer to the next column.
    constexpr void next_column() noexcept
    {
        base_t::next_column();
        trace_matrix_iter += band_dimension - ((current_matrix_iter != seqan3::begin(score_matrix)) ? 1 : 0);
    }

    //!\brief The data container.
    trace_matrix_type trace_matrix{};
    //!\brief The current iterator in the trace matrix.
    typename std::ranges::iterator_t<trace_matrix_type> trace_matrix_iter{};
    //!\brief The full dimension of the band.
    uint_fast32_t band_dimension{};
};
} // namespace seqan3::detail
