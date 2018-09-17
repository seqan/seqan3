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
 * \brief Contains the declaration of seqan3::detail::alignment_trace_matrix.
 */

#pragma once

#include <deque>
#include <vector>

#include <seqan3/alignment/matrix/matrix_concept.hpp>
#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/core/add_enum_bitwise_operators.hpp>
#include <seqan3/core/metafunction/range.hpp>

namespace seqan3::detail
{

/*!\brief The possible directions a trace can have. The values can be combined by the logical `|`-operator.
 * \ingroup alignment_matrix
 * \sa seqan3::add_enum_bitwise_operators <seqan3::detail::trace_directions> enables combining enum values.
 * \sa seqan3::detail::alignment_trace_matrix implementations use this enum as matrix entry type.
 */
enum struct trace_directions : uint8_t
{
    //!\brief No trace
    none      = 0b0000,
    //!\brief Trace comes from the diagonal entry.
    diagonal  = 0b0001,
    //!\brief Trace comes from the above entry.
    up        = 0b0010,
    //!\brief Trace comes from the left entry.
    left      = 0b0100
};

} // namespace seqan3::detail

namespace seqan3
{
//!\brief Enable bitwise operators for enum seqan3::detail::trace_directions.
//!\ingroup alignment_matrix
template <>
constexpr bool add_enum_bitwise_operators<seqan3::detail::trace_directions> = true;
} // namespace seqan3

namespace seqan3::detail
{

//!\brief The declaration of alignment_trace_matrix. Each definition of this
//!       declaration must satisfy seqan3::detail::matrix_concept.
//!\ingroup alignment_matrix
//!\implements seqan3::detail::matrix_concept
template <typename ...>
struct alignment_trace_matrix;

/*!\brief Represents the begin/end of the pairwise alignment in the respective sequences.
 * This class can for example be used to represent the coordinate where the best alignment score is located.
 */
struct alignment_coordinate
{
    //!\brief The position in the first sequence.
    size_t seq1_pos;
    //!\brief The position in the second sequence.
    size_t seq2_pos;
};

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
    using signed_size_t = std::make_signed_t<size_t>;

    constexpr auto N = trace_directions::none;
    constexpr auto D = trace_directions::diagonal;
    constexpr auto L = trace_directions::left;
    constexpr auto U = trace_directions::up;
    signed_size_t row = end_coordinate.seq2_pos + 1;
    signed_size_t col = end_coordinate.seq1_pos + 1;

    assert(row < matrix.rows());
    assert(col < matrix.cols());

    while (true)
    {
        trace_directions dir = matrix.at(row, col);
        if ((dir & L) == L)
        {
            col = std::max<signed_size_t>(col - 1, 0);
        }
        else if ((dir & U) == U)
        {
            row = std::max<signed_size_t>(row - 1, 0);
        }
        else if ((dir & D) == D)
        {
            row = std::max<signed_size_t>(row - 1, 0);
            col = std::max<signed_size_t>(col - 1, 0);
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

    return {std::max<signed_size_t>(col - 1, 0), std::max<signed_size_t>(row - 1, 0)};
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
    using signed_size_t = std::make_signed_t<size_t>;

    constexpr auto N = trace_directions::none;
    constexpr auto D = trace_directions::diagonal;
    constexpr auto L = trace_directions::left;
    constexpr auto U = trace_directions::up;
    signed_size_t col = end_coordinate.seq1_pos + 1;
    signed_size_t row = end_coordinate.seq2_pos + 1;

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
            col = std::max<signed_size_t>(col - 1, 0);
            gapped_database.push_front(database[col]);
            gapped_query.push_front(gap::GAP);
        }
        else if ((dir & U) == U)
        {
            row = std::max<signed_size_t>(row - 1, 0);
            gapped_database.push_front(gap::GAP);
            gapped_query.push_front(query[row]);
        }
        else if ((dir & D) == D)
        {
            row = std::max<signed_size_t>(row - 1, 0);
            col = std::max<signed_size_t>(col - 1, 0);
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

/*!\brief A trace matrix represented in a one-dimensional std::vector
 * \ingroup alignment_matrix
 *
 * \details
 *
 * This data structure stores the matrix in a flat way using the
 * std::vector<#entry_type> data structure where each row is stored
 * continuously.
 *
 * ## Example
 *
 * \snippet test/snippet/alignment/matrix/alignment_trace_matrix_vector.cpp code
 *
 * ### Output
 * \include test/snippet/alignment/matrix/alignment_trace_matrix_vector.out
 */
template <>
struct alignment_trace_matrix<std::vector<trace_directions>>
    : public row_wise_matrix<trace_directions>
{
    using row_wise_matrix<trace_directions>::row_wise_matrix;
};

/*!\brief A trace matrix that uses an underlying seqan3::detail::alignment_score_matrix
 * \ingroup alignment_matrix
 * \tparam database_type            The type of the database sequence.
 * \tparam query_type               The type of the query sequence.
 * \tparam align_config_type        The type of the alignment config.
 * \tparam ...score_matrix_params_t The template parameters of seqan3::detail::alignment_score_matrix
 *
 * \details
 *
 * This data structure uses directly the score matrix to infer the trace matrix
 * and works for any seqan3::detail::alignment_score_matrix.
 *
 * \todo TODO: Is currently only able to handle the edit distance.
 *
 * ## Example
 *
 * \snippet test/snippet/alignment/matrix/alignment_trace_matrix.cpp code
 *
 * ### Output
 * \include test/snippet/alignment/matrix/alignment_trace_matrix.out
 */
template <typename database_type, typename query_type, typename align_config_type, typename ...score_matrix_params_t>
//!\cond
    requires matrix_concept<alignment_score_matrix<score_matrix_params_t...>> &&
             std::Integral<typename alignment_score_matrix<score_matrix_params_t...>::entry_type>
//!\endcond
struct alignment_trace_matrix<database_type, query_type, align_config_type, alignment_score_matrix<score_matrix_params_t...>>
    : public alignment_score_matrix<score_matrix_params_t...>
{
    //!\brief The type of the score matrix.
    using score_matrix_type = alignment_score_matrix<score_matrix_params_t...>;

    //!\brief The type of an entry in the score matrix.
    using score_type = typename score_matrix_type::entry_type;

    //!\copydoc seqan3::detail::matrix_concept::entry_type
    using entry_type = trace_directions;

    /*!\name Constructors, destructor and assignment
     * The copy-constructor, move-constructor, copy-assignment, move-assignment,
     * and destructor are implicitly defined.
     * \{
     */
     alignment_trace_matrix() = default;
     alignment_trace_matrix(alignment_trace_matrix const &) = default;
     alignment_trace_matrix(alignment_trace_matrix &&) = default;
     alignment_trace_matrix & operator=(alignment_trace_matrix const &) = default;
     alignment_trace_matrix & operator=(alignment_trace_matrix &&) = default;

    /*!\brief Construct the trace matrix by using a score_matrix.
     * \param database     The database sequence.
     * \param query        The query sequence.
     * \param config       The alignment config.
     * \param score_matrix The score matrix.
     */
    alignment_trace_matrix(database_type database, query_type query, align_config_type config, score_matrix_type score_matrix)
        : score_matrix_type(std::move(score_matrix)),
          _database{std::forward<database_type>(database)},
          _query{std::forward<query_type>(query)},
          _config{std::forward<align_config_type>(config)}
    {}
    //!\}

    //!\copydoc seqan3::detail::matrix_concept::rows
    using score_matrix_type::rows;
    //!\copydoc seqan3::detail::matrix_concept::cols
    using score_matrix_type::cols;

    //!\brief The trace directions of the matrix at position (*row*, *col*).
    entry_type at(size_t const row, size_t const col) const noexcept
    {
        entry_type direction{};

        if (is_trace_diagonal(row, col))
            direction |= entry_type::diagonal;

        if (is_trace_up(row, col))
            direction |= entry_type::up;

        if (is_trace_left(row, col))
            direction |= entry_type::left;

        return direction;
    }

    //!\brief Access to the score_matrix.
    score_matrix_type const & score_matrix() const noexcept
    {
        return *this;
    }

private:

    //!\brief Does the trace come from the above entry?
    bool is_trace_up(size_t const row, size_t const col) const noexcept
    {
        // TODO: use the alignment_config to calculate the score
        score_type gap = 1;

        score_type curr = score_matrix().at(row, col);
        score_type up = row == 0 ? col : score_matrix().at(row - 1, col);
        return curr == up + gap;
    }

    //!\brief Does the trace come from the left entry?
    bool is_trace_left(size_t const row, size_t const col) const noexcept
    {
        // TODO: use the alignment_config to calculate the score
        score_type gap = 1;

        score_type curr = score_matrix().at(row, col);
        score_type left = col == 0 ? row : score_matrix().at(row, col - 1);
        return curr == left + gap;
    }

    //!\brief Does the trace come from the diagonal entry?
    bool is_trace_diagonal(size_t const row, size_t const col) const noexcept
    {
        // TODO: use the alignment_config to calculate the score
        score_type match = 0;
        score_type mismatch = 1;

        score_type curr = score_matrix().at(row, col);
        if (col == 0 || row == 0)
            return false;

        score_type diag = score_matrix().at(row - 1, col - 1);
        bool is_match = _query[row - 1] == _database[col - 1];

        return (is_match && curr == diag + match) ||
              (!is_match && curr == diag + mismatch);
    }

    //!\brief The database sequence.
    database_type _database;
    //!\brief The query sequence.
    query_type _query;
    //!\brief The alignment config.
    align_config_type _config;
};

/*!\name Type deduction guides
 * \relates seqan3::detail::alignment_trace_matrix
 * \{
 */
alignment_trace_matrix(std::vector<trace_directions>, size_t rows, size_t cols)
    -> alignment_trace_matrix<std::vector<trace_directions>>;

template <typename database_t, typename query_t, typename align_config_t, typename alignment_t, typename ...options_t>
alignment_trace_matrix(database_t && database, query_t && query, align_config_t && config, alignment_score_matrix<alignment_t, options_t...>)
    -> alignment_trace_matrix<database_t, query_t, align_config_t, alignment_score_matrix<alignment_t, options_t...>>;
//!\}

} // namespace seqan3::detail
