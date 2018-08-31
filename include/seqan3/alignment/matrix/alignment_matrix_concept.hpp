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
 * \brief Contains seqan3::alignment_matrix_concept.
 */

#pragma once

#include <cstddef>

#include <seqan3/std/concept/comparison.hpp>

namespace seqan3
{

/*!\interface seqan3::alignment_matrix_concept <>
 * \brief Defines the requirements of an alignment matrix (e.g. score, trace matrices).
 * \tparam matrix_t The type the concept check is performed on (the putative alignment matrix).
 * \ingroup alignment_matrix
 */
/*!\name Requirements for seqan3::alignment_matrix_concept
* \brief You can expect these members on all types that implement seqan3::alignment_matrix_concept.
* \memberof seqan3::alignment_matrix_concept
* \{
*/
//!\cond
template <typename matrix_t>
concept bool alignment_matrix_concept = requires(matrix_t m)
{
//!\endcond
    /*!\typedef typedef IMPLEMENTATION_DEFINED sequence_type;
     * \brief The type of the database and query sequence.
     * \memberof seqan3::alignment_matrix_concept
     */
    typename matrix_t::sequence_type;

    /*!\typedef typedef IMPLEMENTATION_DEFINED entry_type;
     * \brief The type of an entry in the matrix.
     * \memberof seqan3::alignment_matrix_concept
     */
    typename matrix_t::entry_type;

    /*!\fn sequence_type const & database() const noexcept;
     * \brief The database sequence (sequence at the top of the matrix).
     * \memberof seqan3::alignment_matrix_concept
     */
    { m.database() } -> typename matrix_t::sequence_type const &;

    /*!\fn sequence_type const & query() const noexcept;
     * \brief The query sequence (sequence to the left of the matrix).
     * \memberof seqan3::alignment_matrix_concept
     */
    { m.query() } -> typename matrix_t::sequence_type const &;

    /*!\fn std::size_t cols() const noexcept;
     * \brief The number of columns in the matrix.
     * \memberof seqan3::alignment_matrix_concept
     */
    { m.cols() } -> std::size_t;

    /*!\fn std::size_t rows() const noexcept;
     * \brief The number of rows in the matrix.
     * \memberof seqan3::alignment_matrix_concept
     */
    { m.rows() } -> std::size_t;

    /*!\fn entry_type at(unsigned row, unsigned col) const noexcept;
     * \brief The entry of the matrix at position (\a row, \a col), e.g. `matrix[row][col]`.
     * \memberof seqan3::alignment_matrix_concept
     */
    { m.at(0u, 0u) } -> typename matrix_t::entry_type;
//!\cond
};
//!\endcond
//!\}

/*!\brief `true`, if two alignment matrices are equal.
 * \relates alignment_matrix_concept
 * \tparam    matrix1_t The type of the left hand side alignment matrix.
 * \tparam    matrix2_t The type of the right hand side alignment matrix.
 * \param[in] lhs       Compare the left hand side matrix...
 * \param[in] rhs       ... with the right hand side matrix.
 */
template <alignment_matrix_concept matrix1_t, alignment_matrix_concept matrix2_t>
//!\cond
    requires equality_comparable_with_concept<
                typename matrix1_t::entry_type,
                typename matrix2_t::entry_type
             >
//!\endcond
inline bool operator==(matrix1_t const & lhs, matrix2_t const & rhs)
{
    if (lhs.rows() != rhs.rows())
        return false;

    if (lhs.cols() != rhs.cols())
        return false;

    for (unsigned row = 0u; row < lhs.rows(); ++row)
        for (unsigned col = 0u; col < lhs.cols(); ++col)
            if (!(lhs.at(row, col) == rhs.at(row, col)))
                return false;

    return true;
}

/*!\brief `true`, if two alignment matrices are unequal.
 * \relates alignment_matrix_concept
 * \tparam    matrix1_t The type of the left hand side alignment matrix.
 * \tparam    matrix2_t The type of the right hand side alignment matrix.
 * \param[in] lhs       Compare the left hand side matrix...
 * \param[in] rhs       ... with the right hand side matrix.
 */
template <alignment_matrix_concept matrix1_t, alignment_matrix_concept matrix2_t>
//!\cond
    requires equality_comparable_with_concept<
                typename matrix1_t::entry_type,
                typename matrix2_t::entry_type
             >
//!\endcond
inline bool operator!=(matrix1_t const & lhs, matrix2_t const & rhs)
{
    return !(lhs == rhs);
}

} // namespace seqan3
