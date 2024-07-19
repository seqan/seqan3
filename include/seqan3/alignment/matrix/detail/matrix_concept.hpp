// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::detail::matrix.
 */

#pragma once

#include <concepts>
#include <cstddef>
#include <limits>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>

namespace seqan3::detail
{

//!\brief A special score which represents infinity.
//!\ingroup alignment_matrix
template <typename score_type>
constexpr score_type matrix_inf = std::numeric_limits<score_type>::max();

/*!\interface seqan3::detail::matrix <>
 * \brief Defines the requirements of a matrix (e.g. score matrices, trace matrices).
 * \tparam matrix_t The type the concept check is performed on (the putative matrix).
 * \ingroup alignment_matrix
 */
//!\cond
template <typename matrix_t>
concept matrix = requires (std::remove_cvref_t<matrix_t> m) {
    typename std::remove_cvref_t<matrix_t>::value_type;

    typename std::remove_cvref_t<matrix_t>::reference;

    typename std::remove_cvref_t<matrix_t>::size_type;

    { m.cols() } -> std::same_as<typename std::remove_cvref_t<matrix_t>::size_type>;

    { m.rows() } -> std::same_as<typename std::remove_cvref_t<matrix_t>::size_type>;

    { m.at(matrix_coordinate{}) } -> std::same_as<typename std::remove_cvref_t<matrix_t>::reference>;
};
//!\endcond

// Workaround for https://github.com/doxygen/doxygen/issues/9379
#if SEQAN3_DOXYGEN_ONLY(1) 0
template <typename matrix_t>
class matrix
{};
#endif

/*!\name Requirements for seqan3::detail::matrix
 * \brief You can expect these members on all types that implement seqan3::detail::matrix.
 * \relates seqan3::detail::matrix
 * \{
 */

/*!\typedef typedef IMPLEMENTATION_DEFINED value_type;
 * \brief The type of an entry in the matrix.
 */
/*!\typedef typedef IMPLEMENTATION_DEFINED reference;
 * \brief The type of a reference to an entry in the matrix.
 */
/*!\typedef typedef IMPLEMENTATION_DEFINED size_type;
 * \brief The size type of the matrix.
 */
/*!\fn size_type cols() const noexcept;
 * \brief The number of columns in the matrix.
 */
/*!\fn size_type rows() const noexcept;
 * \brief The number of rows in the matrix.
 */
/*!\fn reference at(matrix_coordinate coordinate) noexcept;
 * \brief A reference to the entry of the matrix at the given coordinate.
 */

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
template <matrix matrix1_t, matrix matrix2_t>
    requires std::equality_comparable_with<typename matrix1_t::reference, typename matrix2_t::reference>
inline bool operator==(matrix1_t const & lhs, matrix2_t const & rhs) noexcept
{
    if (lhs.rows() != rhs.rows())
        return false;

    if (lhs.cols() != rhs.cols())
        return false;

    for (size_t row = 0u; row < lhs.rows(); ++row)
        for (size_t col = 0u; col < lhs.cols(); ++col)
            if (matrix_coordinate co{row_index_type{row}, column_index_type{col}}; lhs.at(co) != rhs.at(co))
                return false;

    return true;
}

/*!\brief Whether two alignment matrices are equal.
 * \tparam    matrix1_t The type of the left hand side matrix.
 * \tparam    matrix2_t The type of the right hand side matrix.
 * \param[in] lhs       Compare the left hand side matrix
 * \param[in] rhs       with the right hand side matrix.
 */
template <matrix matrix1_t, matrix matrix2_t>
    requires std::equality_comparable_with<typename matrix1_t::reference, typename matrix2_t::reference>
inline bool operator!=(matrix1_t const & lhs, matrix2_t const & rhs) noexcept
{
    return !(lhs == rhs);
}
//!\}

} // namespace seqan3::detail
