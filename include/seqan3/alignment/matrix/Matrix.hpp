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
 * \brief Contains seqan3::detail::Matrix.
 */

#pragma once

#include <cstddef>

#include <seqan3/std/concepts>

namespace seqan3::detail
{

/*!\interface seqan3::detail::Matrix <>
 * \brief Defines the requirements of an matrix (e.g. score matrices, trace matrices).
 * \tparam matrix_t The type the concept check is performed on (the putative matrix).
 * \ingroup alignment_matrix
 */
/*!\name Requirements for seqan3::detail::Matrix
* \brief You can expect these members on all types that implement seqan3::detail::Matrix.
* \memberof seqan3::detail::Matrix
* \{
*/
//!\cond
template <typename matrix_t>
concept Matrix = requires(matrix_t m)
{
//!\endcond
    /*!\typedef typedef IMPLEMENTATION_DEFINED entry_type;
     * \brief The type of an entry in the matrix.
     * \memberof seqan3::detail::Matrix
     */
    typename matrix_t::entry_type;

    /*!\fn size_t cols() const noexcept;
     * \brief The number of columns in the matrix.
     * \memberof seqan3::detail::Matrix
     */
    { m.cols() } -> size_t;

    /*!\fn size_t rows() const noexcept;
     * \brief The number of rows in the matrix.
     * \memberof seqan3::detail::Matrix
     */
    { m.rows() } -> size_t;

    /*!\fn entry_type at(size_t row, size_t col) const noexcept;
     * \brief The entry of the matrix at position (\a row, \a col), e.g. `matrix[row][col]`.
     * \memberof seqan3::detail::Matrix
     */
    { m.at(size_t{0u}, size_t{0u}) } -> typename matrix_t::entry_type;
//!\cond
};
//!\endcond
//!\}

/*!\name Comparison operators
 * \{
 */
/*!\brief Whether two alignment matrices are equal.
 * \relates Matrix
 * \tparam    matrix1_t The type of the left hand side matrix.
 * \tparam    matrix2_t The type of the right hand side matrix.
 * \param[in] lhs       Compare the left hand side matrix
 * \param[in] rhs       with the right hand side matrix.
 */
template <Matrix matrix1_t, Matrix matrix2_t>
//!\cond
    requires std::EqualityComparableWith<typename matrix1_t::entry_type, typename matrix2_t::entry_type>
//!\endcond
inline bool operator==(matrix1_t const & lhs, matrix2_t const & rhs) noexcept
{
    if (lhs.rows() != rhs.rows())
        return false;

    if (lhs.cols() != rhs.cols())
        return false;

    for (size_t row = 0u; row < lhs.rows(); ++row)
        for (size_t col = 0u; col < lhs.cols(); ++col)
            if (lhs.at(row, col) != rhs.at(row, col))
                return false;

    return true;
}

/*!\brief Whether two alignment matrices are equal.
 * \relates Matrix
 * \tparam    matrix1_t The type of the left hand side matrix.
 * \tparam    matrix2_t The type of the right hand side matrix.
 * \param[in] lhs       Compare the left hand side matrix
 * \param[in] rhs       with the right hand side matrix.
 */
template <Matrix matrix1_t, Matrix matrix2_t>
//!\cond
    requires std::EqualityComparableWith<typename matrix1_t::entry_type, typename matrix2_t::entry_type>
//!\endcond
inline bool operator!=(matrix1_t const & lhs, matrix2_t const & rhs) noexcept
{
    return !(lhs == rhs);
}
//!\}

} // namespace seqan3
