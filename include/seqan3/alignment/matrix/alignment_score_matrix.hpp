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
 * \brief Contains the declaration of seqan3::alignment_score_matrix.
 */

#pragma once

#include <vector>
#include <seqan3/core/platform.hpp>

namespace seqan3
{

//!\brief The declaration of alignment_score_matrix. Each definition of this
//!       declaration must satisfy seqan3::alignment_matrix_concept.
//!\ingroup alignment_matrix
template <typename ...>
struct alignment_score_matrix;

/*!\brief A score matrix represented in a one-dimensional std::vector
 * \ingroup alignment_matrix
 * \implements seqan3::alignment_matrix_concept
 * \tparam score_t    \copydoc seqan3::alignment_score_matrix<std::vector<score_t>,sequence_t>::entry_type
 * \tparam sequence_t \copydoc seqan3::alignment_matrix_concept::sequence_type
 *
 *
 * This data structure stores the matrix in a flat way using the
 * std::vector<#entry_type> data structure where each row is stored
 * continuously.
 *
 * ## Example
 *
 * \snippet test/example/alignment/matrix/alignment_score_matrix.cpp code
 *
 * ### Output
 * \include test/example/alignment/matrix/alignment_score_matrix.out
 */
template <typename score_t, typename sequence_t>
//!\cond
    requires std::is_integral_v<score_t>
//!\endcond
struct alignment_score_matrix<std::vector<score_t>, sequence_t>
{
    //!\brief The type of the score.
    using entry_type = score_t;

    //!\copydoc seqan3::alignment_matrix_concept::sequence_type
    using sequence_type = sequence_t;

    /*!\name Constructors, destructor and assignment
     * The copy-constructor, move-constructor, copy-assignment, move-assignment,
     * and destructor are implicitly defined.
     * \{
     */
    /*!\brief Construct the score matrix out of the \a scores, the \a database,
     *        and the \a query.
     * \param     scores   The score values as flat #std::vector<#entry_type>.
     * \param[in] database The #database sequence.
     * \param[in] query    The #query sequence.
     *
     * \attention Make sure that the #database and the #query outlives the matrix.
     */
    alignment_score_matrix
    (
        std::vector<entry_type> scores,
        sequence_type const & database,
        sequence_type const & query
    )
        : _scores{std::move(scores)}, _database{database}, _query{query}
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
        return _scores[row * cols() + col];
    }

private:
    //!\brief The matrix as a one-dimensional vector of scores
    //!       (each row is continuously stored).
    std::vector<entry_type> _scores;
    //!\copydoc seqan3::alignment_matrix_concept::database
    sequence_type const & _database;
    //!\copydoc seqan3::alignment_matrix_concept::query
    sequence_type const & _query;
};

/*!\name Type deduction guides
 * \relates seqan3::alignment_score_matrix
 * \{
 */
template <typename score_t, typename sequence_t>
alignment_score_matrix(std::vector<score_t>, sequence_t const &, sequence_t const &)
    -> alignment_score_matrix<std::vector<score_t>, sequence_t>;
//!\}

} // namespace seqan3
