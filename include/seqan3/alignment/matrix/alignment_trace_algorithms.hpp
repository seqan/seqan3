// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides algorithms that operate on seqan3::detail::alignment_trace_matrix.
 */

#pragma once

#include <deque>
#include <vector>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/matrix/matrix_concept.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/views/slice.hpp>

namespace seqan3::detail
{

/*!\brief Compute the front coordinate.
 * \ingroup alignment_matrix
 * \tparam    trace_matrix_t The type of the trace matrix.
 * \param[in] matrix          The trace matrix.
 * \param[in] back_coordinate Where the trace in the matrix ends.
 * \returns Returns the front coordinate.
 */
 template <typename trace_matrix_t>
 //!\cond
     requires matrix<remove_cvref_t<trace_matrix_t>> &&
              std::same_as<typename remove_cvref_t<trace_matrix_t>::value_type, trace_directions>
 //!\endcond
inline alignment_coordinate alignment_front_coordinate(trace_matrix_t && matrix,
                                                       alignment_coordinate const back_coordinate)
{
    constexpr auto D = trace_directions::diagonal;
    constexpr auto L = trace_directions::left;
    constexpr auto U = trace_directions::up;

    matrix_coordinate coordinate{row_index_type{back_coordinate.second}, column_index_type{back_coordinate.first}};

    assert(coordinate.row < matrix.rows());
    assert(coordinate.col < matrix.cols());

    while (true)
    {
        trace_directions dir = matrix.at(coordinate);
        if ((dir & L) == L)
        {
            coordinate.col = std::max<size_t>(coordinate.col, 1) - 1;
        }
        else if ((dir & U) == U)
        {
            coordinate.row = std::max<size_t>(coordinate.row, 1) - 1;
        }
        else if ((dir & D) == D)
        {
            coordinate.row = std::max<size_t>(coordinate.row, 1) - 1;
            coordinate.col = std::max<size_t>(coordinate.col, 1) - 1;
        }
        else
        {
#ifndef NDEBUG
            if (!(coordinate.row == 0 || coordinate.col == 0))
                throw std::logic_error{"Unknown seqan3::trace_direction in an inner cell of the trace matrix."};
#endif
            break;
        }
    }

    return {column_index_type{coordinate.col}, row_index_type{coordinate.row}};
}

/*!\brief Compute the trace from a trace matrix.
 * \ingroup alignment_matrix
 * \tparam    alignment_t                The type of the returned alignment.
 * \tparam    database_t                 The type of the database sequence.
 * \tparam    query_t                    The type of the query sequence.
 * \tparam    trace_matrix_t             The type of the trace matrix.
 * \param[in] database                   The database sequence.
 * \param[in] query                      The query sequence.
 * \param[in] matrix                     The trace matrix.
 * \param[in] back_coordinate            Where the trace in the matrix ends.
 * \param[in] front_coordinate           Where the trace in the matrix starts.
 * \returns Returns a seqan3::aligned_sequence.
 */
template <
    tuple_like alignment_t,
    typename database_t,
    typename query_t,
    typename trace_matrix_t>
//!\cond
    requires matrix<remove_cvref_t<trace_matrix_t>> &&
             std::same_as<typename remove_cvref_t<trace_matrix_t>::value_type, trace_directions> &&
             detail::all_satisfy_aligned_seq<detail::tuple_type_list_t<alignment_t>>
//!\endcond
inline alignment_t alignment_trace(database_t && database,
                                   query_t && query,
                                   trace_matrix_t && matrix,
                                   alignment_coordinate const back_coordinate,
                                   alignment_coordinate const front_coordinate)
{
    constexpr auto N = trace_directions::none;
    constexpr auto D = trace_directions::diagonal;
    constexpr auto L = trace_directions::left;
    constexpr auto U = trace_directions::up;

    matrix_coordinate coordinate{row_index_type{back_coordinate.second}, column_index_type{back_coordinate.first}};

    assert(coordinate.row <= query.size());
    assert(coordinate.col <= database.size());
    assert(coordinate.row < matrix.rows());
    assert(coordinate.col < matrix.cols());

    alignment_t aligned_seq{};
    assign_unaligned(std::get<0>(aligned_seq), views::slice(database, front_coordinate.first, coordinate.col));
    assign_unaligned(std::get<1>(aligned_seq), views::slice(query, front_coordinate.second, coordinate.row));
    auto end_aligned_db = std::ranges::cend(std::get<0>(aligned_seq));
    auto end_aligned_qy = std::ranges::cend(std::get<1>(aligned_seq));

    if (matrix.at({row_index_type{0u}, column_index_type{0u}}) != N)
        throw std::logic_error{"End trace must be NONE"};

    while (true)
    {
        trace_directions dir = matrix.at(coordinate);
        if ((dir & L) == L)
        {
            coordinate.col = std::max<size_t>(coordinate.col, 1) - 1;
            --end_aligned_db;
            end_aligned_qy = insert_gap(std::get<1>(aligned_seq), end_aligned_qy);
        }
        else if ((dir & U) == U)
        {
            coordinate.row = std::max<size_t>(coordinate.row, 1) - 1;
            end_aligned_db = insert_gap(std::get<0>(aligned_seq), end_aligned_db);
            --end_aligned_qy;
        }
        else if ((dir & D) == D)
        {
            coordinate.row = std::max<size_t>(coordinate.row, 1) - 1;
            coordinate.col = std::max<size_t>(coordinate.col, 1) - 1;
            --end_aligned_db;
            --end_aligned_qy;
        }
        else
        {
#ifndef NDEBUG
            if (!(coordinate.row == 0 || coordinate.col == 0))
                throw std::logic_error{"Unknown seqan3::trace_direction in an inner cell of the trace matrix."};
#endif
            break;
        }
    }

    return aligned_seq;
}

} // namespace seqan3::detail
