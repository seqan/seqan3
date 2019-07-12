// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::detail::Matrix.
 */

#pragma once

#include <cstddef>
#include <limits>

#include <seqan3/std/concepts>

namespace seqan3::detail
{

//!\brief A special score which represents infinity.
//!\ingroup alignment_matrix
template <typename score_type>
constexpr score_type matrix_inf = std::numeric_limits<score_type>::max();

/*!\interface seqan3::detail::Matrix <>
 * \brief Defines the requirements of a matrix (e.g. score matrices, trace matrices).
 * \tparam matrix_t The type the concept check is performed on (the putative matrix).
 * \ingroup alignment_matrix
 */
/*!\name Requirements for seqan3::detail::Matrix
 * \brief You can expect these members on all types that implement seqan3::detail::Matrix.
 * \relates seqan3::detail::Matrix
 * \{
 */
//!\cond
template <typename matrix_t>
SEQAN3_CONCEPT Matrix = requires(matrix_t m)
{
//!\endcond

    /*!\typedef typedef IMPLEMENTATION_DEFINED entry_type;
     * \brief The type of an entry in the matrix.
     */
    typename std::remove_reference_t<matrix_t>::entry_type;

    /*!\fn size_t cols() const noexcept;
     * \brief The number of columns in the matrix.
     */
    { m.cols() } -> size_t;

    /*!\fn size_t rows() const noexcept;
     * \brief The number of rows in the matrix.
     */
    { m.rows() } -> size_t;

    /*!\fn entry_type at(size_t row, size_t col) const noexcept;
     * \brief The entry of the matrix at position (\a row, \a col), e.g. `matrix[row][col]`.
     */
    { m.at(size_t{0u}, size_t{0u}) } -> typename std::remove_reference_t<matrix_t>::entry_type;
//!\cond
};
//!\endcond
//!\}

/*!\name Comparison operators
 * \ingroup alignment_matrix
 * \{
 */

/*!\brief Whether two alignment matrices are equal.
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

} // namespace seqan3::detail
