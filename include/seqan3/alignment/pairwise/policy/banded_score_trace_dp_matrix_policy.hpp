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
#include <range/v3/view/zip.hpp>

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/pairwise/policy/banded_score_dp_matrix_policy.hpp>
#include <seqan3/alignment/pairwise/policy/unbanded_score_trace_dp_matrix_policy.hpp>
#include <seqan3/std/span>

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
    using base_t::band_size;

    // Import member functions from base class
    using base_t::current_band_size;
    using base_t::second_range_begin_offset;
    using base_t::band_touches_last_row;
    using base_t::trim_sequences;
    using base_t::map_banded_coordinate_to_range_position;

    /*!\name Member types
     * \{
     */
    //!\brief The underlying cell type of the scoring matrix.
    using cell_type = typename score_allocator_t::value_type;
    //!\brief The underlying cell type of the trace matrix.
    using trace_type = typename trace_allocator_t::value_type;
    //!\brief The type of the score matrix.
    using score_matrix_type = std::vector<cell_type, score_allocator_t>;
    //!\brief The type of the trace matrix.
    using trace_matrix_type = std::vector<trace_type, trace_allocator_t>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr banded_score_trace_dp_matrix_policy() = default;
    constexpr banded_score_trace_dp_matrix_policy(banded_score_trace_dp_matrix_policy const &) = default;
    constexpr banded_score_trace_dp_matrix_policy(banded_score_trace_dp_matrix_policy &&) = default;
    constexpr banded_score_trace_dp_matrix_policy & operator=(banded_score_trace_dp_matrix_policy const &) = default;
    constexpr banded_score_trace_dp_matrix_policy & operator=(banded_score_trace_dp_matrix_policy &&) = default;
    ~banded_score_trace_dp_matrix_policy() = default;
    //!\}

    /*!\brief Allocates the memory for the dynamic programming matrix given the two sequences.
     * \tparam first_range_t   The type of the first sequence (or packed sequences).
     * \tparam second_range_t  The type of the second sequence (or packed sequences).
     * \tparam band_t          The type of the band object.
     * \param[in] first_range  The first sequence (or packed sequences).
     * \param[in] second_range The second sequence (or packed sequences).
     * \param[in] band         The band object.
     */
    template <typename first_range_t, typename second_range_t, typename band_t>
    constexpr void allocate_matrix(first_range_t & first_range, second_range_t & second_range, band_t const & band)
    {
        base_t::allocate_matrix(first_range, second_range, band);

        trace_matrix.resize(band_size * dimension_first_range);
        trace_matrix_iter = std::ranges::begin(trace_matrix) + band_column_index;
    }

    //!\brief Returns the current column of the alignment matrix.
    constexpr auto current_column() noexcept
    {
        auto span = base_t::current_band_size();

        assert(span > 0u);  // The span must always be greater than 0.

        // The begin coordinate in the current column begins at it - begin(matrix).
        // The end coordinate ends at it - begin(matrix) + current_band_size
        advanceable_alignment_coordinate<advanceable_alignment_coordinate_state::row>
            col_begin{column_index_type{current_column_index},
                      row_index_type{static_cast<size_t>(std::ranges::distance(std::ranges::begin(score_matrix),
                                                                               current_matrix_iter))}};
        advanceable_alignment_coordinate<advanceable_alignment_coordinate_state::row>
            col_end{column_index_type{current_column_index}, row_index_type{col_begin.second + span}};

        // Return zip view over current column and current column shifted by one to access the previous horizontal.
        auto zip_score = std::view::zip(std::span{std::addressof(*current_matrix_iter), span},
                                           std::span{std::addressof(*(current_matrix_iter + 1)), span});
        return std::view::zip(std::move(zip_score),
                                 std::view::iota(col_begin, col_end),
                                 std::span{std::addressof(*trace_matrix_iter), span});
    }

    //!\brief Moves internal matrix pointer to the next column.
    constexpr void go_next_column() noexcept
    {
        base_t::go_next_column();
        // Move trace back pointer to begin of next column
        trace_matrix_iter = std::ranges::begin(trace_matrix) + band_size * current_column_index;
        // As long as we shift the band towards the band_column_index, jump to the correct begin.
        if (current_column_index < band_column_index)
            trace_matrix_iter += band_column_index - current_column_index;
    }

    /*!\brief Parses the traceback starting from the given coordinate.
     * \param end_coordinate The coordinate from where to start the traceback.
     *
     * \returns A tuple containing the begin coordinate and a tuple with all seqan3::detail::gap_segment s for the
     *          first sequence and the second sequence.
     */
    constexpr auto parse_traceback(alignment_coordinate const & end_coordinate) const
    {
        // Store the trace segments.
        std::deque<gap_segment> first_segments{};
        std::deque<gap_segment> second_segments{};

        // Put the iterator to the position where the traceback starts.
        auto direction_iter = std::ranges::begin(trace_matrix);
        std::ranges::advance(direction_iter, end_coordinate.first * band_size +
                                     end_coordinate.second);

        // Parse the trace until interrupt.
        while (*direction_iter != trace_directions::none)
        {
            // parse until end of diagonal run
            while (static_cast<bool>(*direction_iter & trace_directions::diagonal))
            {
                std::ranges::advance(direction_iter, -static_cast<int_fast32_t>(band_size));
            }

            size_t col_pos = std::ranges::distance(std::ranges::begin(trace_matrix), direction_iter) / band_size;
            // parse vertical gap -> record gap in first_segments (will be translated into gap of first sequence)
            if (static_cast<bool>(*direction_iter & trace_directions::up) ||
                static_cast<bool>(*direction_iter & trace_directions::up_open))
            {
                // Get the current column index (note the column based layout)

                gap_segment gap{col_pos, 0u};

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
                size_t pos = std::ranges::distance(std::ranges::begin(trace_matrix), direction_iter) % band_size +
                             static_cast<int_fast32_t>(col_pos - band_column_index);
                gap_segment gap{pos, 0u};

                // Follow gap until open signal is detected.
                while (!static_cast<bool>(*direction_iter & trace_directions::left_open))
                {
                    std::ranges::advance(direction_iter, -static_cast<int_fast32_t>(band_size) + 1);
                    ++gap.size;
                }
                // explicitly follow opening gap
                std::ranges::advance(direction_iter, -static_cast<int_fast32_t>(band_size) + 1);
                ++gap.size;
                second_segments.push_front(std::move(gap));
            }
        }

        // Get begin coordinate.
        auto c = column_index_type{
                static_cast<uint_fast32_t>(std::ranges::distance(std::ranges::begin(trace_matrix), direction_iter) /
                                           band_size)};
        auto r = row_index_type{
                static_cast<uint_fast32_t>(std::ranges::distance(std::ranges::begin(trace_matrix), direction_iter) %
                                           band_size)};

        // Validate correct coordinates.
        auto begin_coordinate = map_banded_coordinate_to_range_position(
                alignment_coordinate{column_index_type{c}, row_index_type{r}});
        assert(begin_coordinate.first >= 0u);
        assert(begin_coordinate.first <= end_coordinate.first);

        assert(begin_coordinate.second >= 0u);
        assert(begin_coordinate.second <= map_banded_coordinate_to_range_position(end_coordinate).second);

        return std::tuple{begin_coordinate, first_segments, second_segments};
    }

    //!\brief Helper function to print the trace matrix; for debugging only.
    constexpr void print_trace_matrix() const
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

        // First part: moving band right
        for (size_t col = 0; col < band_column_index; ++col)
        {
            auto it = std::ranges::begin(trace_matrix) + (band_size * col) + (band_column_index - col);
            for (size_t row = 0; row <= std::min(dimension_second_range - 1, band_row_index + col); ++row, ++it)
            {
                debug_stream << /*std::ranges::distance(std::ranges::begin(trace_matrix), it)*/printable(*it) << ',';
            }
            debug_stream << '\n';
        }

        // Second part: moving band down.
        for (size_t col = band_column_index; col < dimension_first_range; ++col)
        {
            for (size_t padding = 0; padding < col - band_column_index; ++padding)
                debug_stream << " ,";

            auto it = std::ranges::begin(trace_matrix) + (band_size * col);
            for (size_t row = 0; row < band_size; ++row, ++it)
            {
                // If the band moves out of the matrix do not try to print the characters.
                if (col - band_column_index + row >= dimension_second_range)
                    continue;
                debug_stream << printable(*it) << ',';
            }
            debug_stream << '\n';
        }
    }

    //!\brief The data container.
    trace_matrix_type trace_matrix{};
    //!\brief The current iterator in the trace matrix.
    typename std::ranges::iterator_t<trace_matrix_type> trace_matrix_iter{};
};
} // namespace seqan3::detail
