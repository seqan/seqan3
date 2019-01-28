// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Contains seqan3::detail::matrix_concept.
 */

#pragma once

#include <cstddef>

#include <seqan3/std/concepts>

namespace seqan3::detail
{

/*!\interface seqan3::detail::matrix_concept <>
 * \brief Defines the requirements of an matrix (e.g. score matrices, trace matrices).
 * \tparam matrix_t The type the concept check is performed on (the putative matrix).
 * \ingroup alignment_matrix
 */
/*!\name Requirements for seqan3::detail::matrix_concept
* \brief You can expect these members on all types that implement seqan3::detail::matrix_concept.
* \memberof seqan3::detail::matrix_concept
* \{
*/
//!\cond
template <typename matrix_t>
SEQAN3_CONCEPT matrix_concept = requires(matrix_t m)
{
//!\endcond
    /*!\typedef typedef IMPLEMENTATION_DEFINED entry_type;
     * \brief The type of an entry in the matrix.
     * \memberof seqan3::detail::matrix_concept
     */
    typename matrix_t::entry_type;

    /*!\fn size_t cols() const noexcept;
     * \brief The number of columns in the matrix.
     * \memberof seqan3::detail::matrix_concept
     */
    { m.cols() } -> size_t;

    /*!\fn size_t rows() const noexcept;
     * \brief The number of rows in the matrix.
     * \memberof seqan3::detail::matrix_concept
     */
    { m.rows() } -> size_t;

    /*!\fn entry_type at(size_t row, size_t col) const noexcept;
     * \brief The entry of the matrix at position (\a row, \a col), e.g. `matrix[row][col]`.
     * \memberof seqan3::detail::matrix_concept
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
 * \relates matrix_concept
 * \tparam    matrix1_t The type of the left hand side matrix.
 * \tparam    matrix2_t The type of the right hand side matrix.
 * \param[in] lhs       Compare the left hand side matrix
 * \param[in] rhs       with the right hand side matrix.
 */
template <matrix_concept matrix1_t, matrix_concept matrix2_t>
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
 * \relates matrix_concept
 * \tparam    matrix1_t The type of the left hand side matrix.
 * \tparam    matrix2_t The type of the right hand side matrix.
 * \param[in] lhs       Compare the left hand side matrix
 * \param[in] rhs       with the right hand side matrix.
 */
template <matrix_concept matrix1_t, matrix_concept matrix2_t>
//!\cond
    requires std::EqualityComparableWith<typename matrix1_t::entry_type, typename matrix2_t::entry_type>
//!\endcond
inline bool operator!=(matrix1_t const & lhs, matrix2_t const & rhs) noexcept
{
    return !(lhs == rhs);
}
//!\}

} // namespace seqan3
