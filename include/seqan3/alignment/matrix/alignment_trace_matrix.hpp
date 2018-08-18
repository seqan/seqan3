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
 * \brief Contains the declaration of seqan3::alignment_trace_matrix.
 */

#pragma once

#include <seqan3/alignment/matrix/alignment_matrix_concept.hpp>
#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/core/add_enum_bitwise_operators.hpp>

namespace seqan3
{

/*!\brief The possible directions a trace can have. The values can be combined.
 * \ingroup alignment_matrix
 * \sa seqan3::add_enum_bitwise_operators<trace_matrix_directions> enables combining enum values.
 * \sa seqan3::alignment_trace_matrix implementations use this enum as matrix entry type.
 */
enum struct trace_matrix_directions : uint8_t
{
    //!\brief No trace
    none      = 0,
    //!\brief Trace comes from the diagonal entry.
    diagonal  = 1,
    //!\brief Trace comes from the above entry.
    up        = 2,
    //!\brief Trace comes from the left entry.
    left      = 4
};

//!\brief Enable bitwise operators for enum seqan3::trace_matrix_directions.
//!\ingroup alignment_matrix
template <>
constexpr bool add_enum_bitwise_operators<trace_matrix_directions> = true;

//!\brief The declaration of alignment_trace_matrix. Each definition of this
//!       declaration must satisfy seqan3::alignment_matrix_concept.
//!\ingroup alignment_matrix
template <typename ...>
struct alignment_trace_matrix;

/*!\brief A trace matrix represented in a one-dimensional std::vector
 * \ingroup alignment_matrix
 * \implements seqan3::alignment_matrix_concept
 * \tparam sequence_t \copydoc seqan3::alignment_matrix_concept::sequence_type
 *
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
template <typename sequence_t>
struct alignment_trace_matrix<std::vector<trace_matrix_directions>, sequence_t>
{
    //!\brief The type of the trace.
    using entry_type = trace_matrix_directions;

    //!\copydoc seqan3::alignment_matrix_concept::sequence_type
    using sequence_type = sequence_t;

    /*!\name Constructors, destructor and assignment
     * The copy-constructor, move-constructor, copy-assignment, move-assignment,
     * and destructor are implicitly defined.
     * \{
     */
    /*!\brief Construct the trace matrix out of the \a traces, the \a database,
     *        and the \a query.
     * \param     traces   The trace values as flat std::vector<#entry_type>.
     * \param[in] database The #database sequence.
     * \param[in] query    The #query sequence.
     *
     * \attention Make sure that the #database and the #query outlives the matrix.
     */
    alignment_trace_matrix
    (
        std::vector<entry_type> traces,
        sequence_type const & database,
        sequence_type const & query
    )
        : _traces{std::move(traces)}, _database{database}, _query{query}
    {}
    //!\}

    //!\copydoc seqan3::alignment_matrix_concept::database
    inline sequence_type const & database() const noexcept
    {
        return _database;
    }

    //!\copydoc seqan3::alignment_matrix_concept::query
    inline sequence_type const & query() const noexcept
    {
        return _query;
    }

    //!\copydoc seqan3::alignment_matrix_concept::rows
    inline std::size_t rows() const noexcept
    {
        return _query.size()+1;
    }

    //!\copydoc seqan3::alignment_matrix_concept::cols
    inline std::size_t cols() const noexcept
    {
        return _database.size()+1;
    }

    //!\copydoc seqan3::alignment_matrix_concept::at
    inline entry_type at(unsigned row, unsigned col) const noexcept
    {
        return _traces[row * cols() + col];
    }

private:
    //!\brief The matrix as a one-dimensional vector of traces
    //!       (each row is continuously stored).
    std::vector<entry_type> _traces;
    //!\copydoc seqan3::alignment_matrix_concept::database
    sequence_type const & _database;
    //!\copydoc seqan3::alignment_matrix_concept::query
    sequence_type const & _query;
};

/*!\brief A trace matrix that uses an underlying seqan3::alignment_score_matrix
 * \ingroup alignment_matrix
 * \implements seqan3::alignment_matrix_concept
 * \tparam ...score_matrix_params_t The template parameters of seqan3::alignment_score_matrix
 *
 *
 * This data structure uses directly the score matrix to infer the trace matrix
 * and works for any seqan3::alignment_score_matrix.
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
template <typename ...score_matrix_params_t>
//!\cond
    requires alignment_matrix_concept<alignment_score_matrix<score_matrix_params_t...>> &&
             std::is_integral_v<typename alignment_score_matrix<score_matrix_params_t...>::entry_type>
//!\endcond
struct alignment_trace_matrix<alignment_score_matrix<score_matrix_params_t...>>
    : public alignment_score_matrix<score_matrix_params_t...>
{
    //!\brief The type of the score matrix.
    using score_matrix_type = alignment_score_matrix<score_matrix_params_t...>;

    //!\brief The type of an entry in the score matrix.
    using score_type = typename score_matrix_type::entry_type;

    //!\copydoc seqan3::alignment_matrix_concept::entry_type
    using entry_type = trace_matrix_directions;

    //!\copydoc seqan3::alignment_matrix_concept::sequence_type
    using sequence_type = typename score_matrix_type::sequence_type;

    /*!\name Constructors, destructor and assignment
     * The copy-constructor, move-constructor, copy-assignment, move-assignment,
     * and destructor are implicitly defined.
     * \{
     */
    //!\brief Construct the trace matrix by using a score_matrix.
    //!\param score_matrix The score matrix.
    alignment_trace_matrix(score_matrix_type score_matrix)
        : score_matrix_type(std::move(score_matrix))
    {}

    //!\brief Construct the trace matrix directly via the score matrix constructor.
    //!\tparam ...args_t The type of the arguments needed to construct the #score_matrix_type
    //!\param  args      The constructor arguments of the base class #score_matrix_type
    template <typename ...args_t>
    alignment_trace_matrix(args_t && ...args)
        : score_matrix_type(std::forward<args_t>(args)...)
    {}
    //!\}

    //!\copydoc seqan3::alignment_matrix_concept::rows
    using score_matrix_type::rows;
    //!\copydoc seqan3::alignment_matrix_concept::cols
    using score_matrix_type::cols;
    //!\copydoc seqan3::alignment_matrix_concept::database
    using score_matrix_type::database;
    //!\copydoc seqan3::alignment_matrix_concept::query
    using score_matrix_type::query;

    //!\brief The trace directions of the matrix at position (\a row, \a col).
    inline entry_type at(unsigned row, unsigned col) const noexcept
    {
        entry_type direction{};

        if(is_trace_diagonal(row, col))
            direction |= entry_type::diagonal;

        if(is_trace_up(row, col))
            direction |= entry_type::up;

        if(is_trace_left(row, col))
            direction |= entry_type::left;

        return direction;
    }

    //!\brief Access to the score_matrix.
    inline score_matrix_type const & score_matrix() const noexcept
    {
        return *this;
    }

private:

    //!\brief Does the trace come from the above entry?
    inline bool is_trace_up(unsigned row, unsigned col) const noexcept
    {
        // TODO: use the alignment_config to calculate the score
        score_type gap = 1;

        score_type curr = score_matrix().at(row, col);
        score_type up = row == 0 ? col : score_matrix().at(row-1, col);
        return curr == up + gap;
    }

    //!\brief Does the trace come from the left entry?
    inline bool is_trace_left(unsigned row, unsigned col) const noexcept
    {
        // TODO: use the alignment_config to calculate the score
        score_type gap = 1;

        score_type curr = score_matrix().at(row, col);
        score_type left = col == 0 ? row : score_matrix().at(row, col-1);
        return curr == left + gap;
    }

    //!\brief Does the trace come from the diagonal entry?
    inline bool is_trace_diagonal(unsigned row, unsigned col) const noexcept
    {
        // TODO: use the alignment_config to calculate the score
        score_type match = 0;
        score_type mismatch = 1;

        score_type curr = score_matrix().at(row, col);
        if (col == 0 || row == 0)
            return false;

        score_type diag = score_matrix().at(row-1, col-1);
        bool is_match = query()[row-1] == database()[col-1];

        return (is_match && curr == diag + match) ||
              (!is_match && curr == diag + mismatch);
    }
};

/*!\name Type deduction guides
 * \relates seqan3::alignment_trace_matrix
 * \{
 */
template <typename sequence_t>
alignment_trace_matrix(std::vector<trace_matrix_directions>, sequence_t const &, sequence_t const &)
    -> alignment_trace_matrix<std::vector<trace_matrix_directions>, sequence_t const &>;

template <typename ...args_t>
alignment_trace_matrix(args_t&& ...args)
    -> alignment_trace_matrix<decltype(alignment_score_matrix{std::forward<args_t>(args)...})>;

template <typename alignment_t, typename ...options_t>
alignment_trace_matrix(alignment_score_matrix<alignment_t, options_t...>)
    -> alignment_trace_matrix<alignment_score_matrix<alignment_t, options_t...>>;
//!\}

} // namespace seqan3
