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

#include <deque>

#include <range/v3/view/iota.hpp>
#include <range/v3/view/reverse.hpp>
#include <range/v3/view/zip.hpp>

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alignment/pairwise/policy/unbanded_dp_matrix_policy.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/std/span.hpp>

namespace seqan3::detail
{

/*!\brief Stores information about a contiguous gap.
 * \ingroup alignment_policy
 */
struct gap_segment
{
    //!\brief The position in the sequence where to insert the gap; (inserted before the position).
    size_t position;
    //!\brief The size of the gap.
    size_t size;
};

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
    //!\brief Make current_column_index visible in this class.
    using base_t::current_column_index;

    /*!\name Member types
     * \{
     */
    //!\brief The underlying cell type of the dynamic programming matrix.
    using cell_type  = typename score_allocator_t::value_type;
    using trace_type = typename trace_allocator_t::value_type;
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
    //!\}

    /*!\brief Allocates the memory for the dynamic programming matrix given the two sequences.
     * \tparam first_batch_t   The type of the first sequence (or packed sequences).
     * \tparam second_batch_t  The type of the second sequence (or packed sequences).
     * \param[in] first_batch  The first sequence (or packed sequences).
     * \param[in] second_batch The second sequence (or packed sequences).
     */
    template <typename first_batch_t, typename second_batch_t>
    constexpr void allocate_matrix(first_batch_t && first_batch, second_batch_t && second_batch)
    {
        dimension_first_batch = seqan3::size(first_batch) + 1;
        dimension_second_batch = seqan3::size(second_batch) + 1;
        current_column_index = 0;
        // We use only one column to compute the score.
        score_matrix.resize(dimension_second_batch);
        // We use the full matrix to store the trace direction.
        trace_matrix.resize(dimension_first_batch * dimension_second_batch);
        trace_matrix_iter = trace_matrix.begin();
    }

    //!\brief Returns the current column of the alignment matrix.
    constexpr auto current_column()
    {
        advanceable_alignment_coordinate<size_t, advanceable_alignment_coordinate_state::row>
            col_begin{column_index_type{current_column_index}, row_index_type{0u}};
        advanceable_alignment_coordinate<size_t, advanceable_alignment_coordinate_state::row>
            col_end{column_index_type{current_column_index}, row_index_type{dimension_second_batch}};

        return ranges::view::zip(std::span{score_matrix},
                                 ranges::view::iota(col_begin, col_end),
                                 std::span{&(*trace_matrix_iter), dimension_second_batch});
    }

    //!\brief Moves internal matrix pointer to the next column.
    constexpr void next_column() noexcept
    {
        ++current_column_index;
        std::advance(trace_matrix_iter, dimension_second_batch);
    }

    /*!\brief Parses the traceback starting from the given coordinate.
     * \param end_coordinate The coordinate from where to start the traceback.
     *
     * \returns A tuple containing the begin coordinate and a tuple with all seqan3::detail::gap_segment s for the
     *          first sequence and the second sequence.
     */
    constexpr auto parse_traceback(alignment_coordinate const & end_coordinate)
    {
        // Store the trace segments.
        std::deque<gap_segment> first_segments{};
        std::deque<gap_segment> second_segments{};

        // Put the iterator to the position where the traceback starts.
        auto direction_iter = seqan3::begin(trace_matrix);
        std::advance(direction_iter, end_coordinate.first_seq_pos * dimension_second_batch +
                                     end_coordinate.second_seq_pos);

        // Parse the trace until interrupt.
        while (*direction_iter != trace_directions::none)
        {
            // parse until end of diagonal run
            while (static_cast<bool>(*direction_iter & trace_directions::diagonal))
            {
                std::advance(direction_iter, -dimension_second_batch - 1);
            }

            // parse vertical gap -> record gap in first_segments (will be translated into gap of first sequence)
            if (static_cast<bool>(*direction_iter & trace_directions::up) ||
                static_cast<bool>(*direction_iter & trace_directions::up_open))
            {
                // Get the current column index (note the column based layout)
                size_t pos = std::distance(seqan3::begin(trace_matrix), direction_iter) / dimension_second_batch;
                gap_segment gap{pos, 0u};

                // Follow gap until open signal is detected.
                while (!static_cast<bool>(*direction_iter & trace_directions::up_open))
                {
                    --direction_iter;
                    ++gap.size;
                }
                // explicitly follow opening gap
                --direction_iter;
                ++gap.size;
                // record the gap
                first_segments.push_front(std::move(gap));
                continue;
            }
            // parse horizontal gap -> record gap in second_segments (will be translated into gap of second sequence)
            if (static_cast<bool>(*direction_iter & trace_directions::left) ||
                static_cast<bool>(*direction_iter & trace_directions::left_open))
            {
                // Get the current row index (note the column based layout)
                size_t pos = std::distance(seqan3::begin(trace_matrix), direction_iter) % dimension_second_batch;
                gap_segment gap{pos, 0u};

                // Follow gap until open signal is detected.
                while (!static_cast<bool>(*direction_iter & trace_directions::left_open))
                {
                    std::advance(direction_iter, -dimension_second_batch);
                    ++gap.size;
                }
                // explicitly follow opening gap
                std::advance(direction_iter, -dimension_second_batch);
                ++gap.size;
                second_segments.push_front(std::move(gap));
            }
        }

        // Get begin coordinate.
        auto c = column_index_type{std::distance(seqan3::begin(trace_matrix), direction_iter) / dimension_second_batch};
        auto r = row_index_type{std::distance(seqan3::begin(trace_matrix), direction_iter) % dimension_second_batch};

        return std::tuple{alignment_coordinate{column_index_type{std::move(c)}, row_index_type{std::move(r)}},
                          first_segments,
                          second_segments};
    }

    //!\brief Prints the trace matrix; for debugging only.
    void print_trace_matrix()
    {
        auto printable = [](trace_directions dir)
        {
            std::string seq{};
            if (dir == trace_directions::none)
                seq.append("0");
            if ((dir & trace_directions::diagonal) == trace_directions::diagonal)
                seq.append("\\");
            if ((dir & trace_directions::up) == trace_directions::up)
                seq.append("|");
            if ((dir & trace_directions::up_open) == trace_directions::up_open)
                seq.append("^");
            if ((dir & trace_directions::left) == trace_directions::left)
                seq.append("-");
            if ((dir & trace_directions::left_open) == trace_directions::left_open)
                seq.append("<");
            return seq;
        };

        for (size_t row = 0; row < dimension_second_batch; ++row)
        {
            for (size_t col = 0; col < dimension_first_batch; ++col)
            {
                debug_stream << printable(trace_matrix[col * dimension_second_batch + row]) << " ";
            }
            debug_stream << "\n";
        }
    }

    //!\brief The data container.
    std::vector<trace_type, trace_allocator_t> trace_matrix{};
    //!\brief The current iterator in the trace matrix.
    typename std::vector<trace_type, trace_allocator_t>::iterator trace_matrix_iter{};
};
} // namespace seqan3::detail
