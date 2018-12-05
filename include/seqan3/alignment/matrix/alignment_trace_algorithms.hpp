// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

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
#include <seqan3/core/metafunction/range.hpp>

namespace seqan3::detail
{

/*!\brief Compute the begin coordinate.
 * \ingroup alignment_matrix
 * \tparam    trace_matrix_t The type of the trace matrix.
 * \param[in] matrix         The trace matrix.
 * \param[in] end_coordinate Where the trace in the matrix ends.
 * \returns Returns the begin coordinate.
 */
 template <typename trace_matrix_t>
 //!\cond
     requires matrix_concept<remove_cvref_t<trace_matrix_t>> &&
              std::Same<typename remove_cvref_t<trace_matrix_t>::entry_type, trace_directions>
 //!\endcond
inline alignment_coordinate alignment_begin_coordinate(trace_matrix_t && matrix,
                                                       alignment_coordinate const end_coordinate)
{
    constexpr auto D = trace_directions::diagonal;
    constexpr auto L = trace_directions::left;
    constexpr auto U = trace_directions::up;
    size_t row = end_coordinate.seq2_pos + 1;
    size_t col = end_coordinate.seq1_pos + 1;

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
                throw std::logic_error{"Unkown seqan3::trace_direction in an inner cell of the trace matrix."};
#endif
            break;
        }
    }

    return {std::max<size_t>(col, 1) - 1, std::max<size_t>(row, 1) - 1};
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
 * \param[in] end_coordinate             Where the trace in the matrix ends.
 * \returns Returns a seqan3::aligned_sequence.
 */
template <
    typename database_t,
    typename query_t,
    typename trace_matrix_t,
    typename gapped_database_alphabet_t = gapped<value_type_t<database_t>>,
    typename gapped_query_alphabet_t = gapped<value_type_t<query_t>>>
//!\cond
    requires matrix_concept<remove_cvref_t<trace_matrix_t>> &&
             std::Same<typename remove_cvref_t<trace_matrix_t>::entry_type, trace_directions>
//!\endcond
inline std::pair<std::vector<gapped_database_alphabet_t>, std::vector<gapped_query_alphabet_t>>
alignment_trace(database_t && database,
                query_t && query,
                trace_matrix_t && matrix,
                alignment_coordinate const end_coordinate)
{
    constexpr auto N = trace_directions::none;
    constexpr auto D = trace_directions::diagonal;
    constexpr auto L = trace_directions::left;
    constexpr auto U = trace_directions::up;
    size_t col = end_coordinate.seq1_pos + 1;
    size_t row = end_coordinate.seq2_pos + 1;

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
            gapped_query.push_front(gap::GAP);
        }
        else if ((dir & U) == U)
        {
            row = std::max<size_t>(row, 1) - 1;
            gapped_database.push_front(gap::GAP);
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
                throw std::logic_error{"Unkown seqan3::trace_direction in an inner cell of the trace matrix."};
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
