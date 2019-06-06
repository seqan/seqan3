// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Contains algorithms that operate on seqan3::detail::alignment_trace_matrix.
 */

#pragma once

#include <deque>
#include <vector>

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/core/type_traits/range.hpp>

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
     requires Matrix<remove_cvref_t<trace_matrix_t>> &&
              std::Same<typename remove_cvref_t<trace_matrix_t>::entry_type, trace_directions>
 //!\endcond
inline alignment_coordinate alignment_front_coordinate(trace_matrix_t && matrix,
                                                       alignment_coordinate const back_coordinate)
{
    constexpr auto D = trace_directions::diagonal;
    constexpr auto L = trace_directions::left;
    constexpr auto U = trace_directions::up;
    size_t row = back_coordinate.second + 1;
    size_t col = back_coordinate.first + 1;

    assert(row < matrix.rows());
    assert(col < matrix.cols());

    while (true)
    {
        trace_directions dir = matrix.at(row, col);
        if ((dir & L) == L)
        {
            col = std::max<size_t>(col, 1) - 1;
        }
        else if ((dir & U) == U)
        {
            row = std::max<size_t>(row, 1) - 1;
        }
        else if ((dir & D) == D)
        {
            row = std::max<size_t>(row, 1) - 1;
            col = std::max<size_t>(col, 1) - 1;
        }
        else
        {
#ifndef NDEBUG
            if (!(row == 0 || col == 0))
                throw std::logic_error{"Unknown seqan3::trace_direction in an inner cell of the trace matrix."};
#endif
            break;
        }
    }

    return {column_index_type{col}, row_index_type{row}};
}

/*!\brief Compute the trace from a trace matrix.
 * \ingroup alignment_matrix
 * \tparam    database_t                 The type of the database sequence.
 * \tparam    query_t                    The type of the query sequence.
 * \tparam    trace_matrix_t             The type of the trace matrix.
 * \cond DEV
 * \tparam    gapped_database_alphabet_t The alphabet type of the gapped database sequence.
 * \tparam    gapped_query_alphabet_t    The alphabet type of the gapped query sequence.
 * \endcond
 * \param[in] database                   The database sequence.
 * \param[in] query                      The query sequence.
 * \param[in] matrix                     The trace matrix.
 * \param[in] back_coordinate            Where the trace in the matrix ends.
 * \returns Returns a seqan3::aligned_sequence.
 */
template <
    typename database_t,
    typename query_t,
    typename trace_matrix_t,
    typename gapped_database_alphabet_t = gapped<value_type_t<database_t>>,
    typename gapped_query_alphabet_t = gapped<value_type_t<query_t>>>
//!\cond
    requires Matrix<remove_cvref_t<trace_matrix_t>> &&
             std::Same<typename remove_cvref_t<trace_matrix_t>::entry_type, trace_directions>
//!\endcond
inline std::pair<std::vector<gapped_database_alphabet_t>, std::vector<gapped_query_alphabet_t>>
alignment_trace(database_t && database,
                query_t && query,
                trace_matrix_t && matrix,
                alignment_coordinate const back_coordinate)
{
    constexpr auto N = trace_directions::none;
    constexpr auto D = trace_directions::diagonal;
    constexpr auto L = trace_directions::left;
    constexpr auto U = trace_directions::up;
    size_t col = back_coordinate.first + 1;
    size_t row = back_coordinate.second + 1;

    assert(row <= query.size());
    assert(col <= database.size());
    assert(row < matrix.rows());
    assert(col < matrix.cols());

    std::deque<gapped_database_alphabet_t> gapped_database{};
    std::deque<gapped_query_alphabet_t> gapped_query{};

    if (matrix.at(0, 0) != N)
        throw std::logic_error{"End trace must be NONE"};

    while (true)
    {
        trace_directions dir = matrix.at(row, col);
        if ((dir & L) == L)
        {
            col = std::max<size_t>(col, 1) - 1;
            gapped_database.push_front(database[col]);
            gapped_query.push_front(gap{});
        }
        else if ((dir & U) == U)
        {
            row = std::max<size_t>(row, 1) - 1;
            gapped_database.push_front(gap{});
            gapped_query.push_front(query[row]);
        }
        else if ((dir & D) == D)
        {
            row = std::max<size_t>(row, 1) - 1;
            col = std::max<size_t>(col, 1) - 1;
            gapped_database.push_front(database[col]);
            gapped_query.push_front(query[row]);
        }
        else
        {
#ifndef NDEBUG
            if (!(row == 0 || col == 0))
                throw std::logic_error{"Unknown seqan3::trace_direction in an inner cell of the trace matrix."};
#endif
            break;
        }
    }

    return
    {
        {std::begin(gapped_database), std::end(gapped_database)},
        {std::begin(gapped_query), std::end(gapped_query)}
    };
}

} // namespace seqan3::detail
