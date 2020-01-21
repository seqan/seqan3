// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_coordinate and associated strong types.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{
/*!\brief A strong type for designated initialisation of the column index of a marix.
 * \ingroup alignment_matrix
 * \tparam index_type The type of the index to store; must model seqan3::arithmetic.
 */
template <std::integral index_type>
struct column_index_type : detail::strong_type<index_type, column_index_type<index_type>>
{
    //!!\brief Import base class constructor.
    using detail::strong_type<index_type, column_index_type<index_type>>::strong_type;
};

/*!\name Type deduction guides
 * \relates column_index_type
 * \{
 */
//!\brief Deduces an signed integral type to std::ptrdiff_t.
template <std::signed_integral index_type>
column_index_type(index_type) -> column_index_type<std::ptrdiff_t>;

//!\brief Deduces an unsigned integral type to size_t.
template <std::unsigned_integral index_type>
column_index_type(index_type) -> column_index_type<size_t>;
//!\}

/*!\brief A strong type for designated initialisation of the row index of a marix.
 * \ingroup alignment_matrix
 * \tparam index_type The type of the index to store; must model seqan3::arithmetic.
 */
template <std::integral index_type>
struct row_index_type : detail::strong_type<index_type, row_index_type<index_type>>
{
    //!!\brief Import base class constructor.
    using detail::strong_type<index_type, row_index_type<index_type>>::strong_type;
};

/*!\name Type deduction guides
 * \relates row_index_type
 * \{
 */
//!\brief Deduces an signed integral type to std::ptrdiff_t.
template <std::signed_integral index_type>
row_index_type(index_type) -> row_index_type<std::ptrdiff_t>;

//!\brief Deduces an unsigned integral type to size_t.
template <std::unsigned_integral index_type>
row_index_type(index_type) -> row_index_type<size_t>;
//!\}

//!\brief A representation of a location or offset within a two-dimensional matrix.
//!\ingroup alignment_matrix
template <std::integral index_t>
struct matrix_index
{
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr matrix_index() = default; //!< Defaulted.
    constexpr matrix_index(matrix_index const &) = default; //!< Defaulted.
    constexpr matrix_index(matrix_index &&) = default; //!< Defaulted.
    constexpr matrix_index & operator=(matrix_index const &) = default; //!< Defaulted.
    constexpr matrix_index & operator=(matrix_index &&) = default; //!< Defaulted.
    ~matrix_index() = default; //!< Defaulted.

    /*!\brief Construction from strongly typed row index and column index.
     * \param row_idx The row index to set.
     * \param col_idx The column index to set.
     */
    constexpr matrix_index(row_index_type<index_t> const row_idx,
                           column_index_type<index_t> const col_idx) noexcept
        : row{row_idx.get()}, col{col_idx.get()}
    {}

    /*!\brief Construction from other matrix_index with different integral type.
     * \param[in] other The other matrix_index to construct from.
     */
    template <std::integral other_index_t>
    //!\cond
        requires !std::same_as<other_index_t, index_t>
    //!\endcond
    explicit constexpr matrix_index(matrix_index<other_index_t> other) noexcept
        : row{static_cast<index_t>(other.row)}, col{static_cast<index_t>(other.col)}
    {}
    //!\}

    //!\brief Explicit conversion to the a std::pair.
    template <std::integral first_index_t, std::integral second_index_t>
    constexpr explicit operator std::pair<first_index_t, second_index_t>() const noexcept
    {
        return std::pair{static_cast<first_index_t>(col), static_cast<second_index_t>(row)};
    }

    index_t row{}; //!< The row index.
    index_t col{}; //!< The column index.
};

/*!\name Type deduction guides
 * \relates matrix_index
 * \{
 */
//!\brief Deduces the default index type to std::ptrdiff_t.
matrix_index() -> matrix_index<std::ptrdiff_t>;

//!\brief Deduces the index type from the common type of both index types.
template <std::integral row_index_t, std::integral col_index_t>
//!\cond
    requires std::common_with<row_index_t, col_index_t>
//!\endcond
matrix_index(row_index_type<row_index_t>, column_index_type<col_index_t>) ->
    matrix_index<std::common_type_t<row_index_t, col_index_t>>;

//!\}

//!\brief A coordinate type to access an element inside of a two-dimensional matrix.
//!\ingroup alignment_matrix
using matrix_coordinate = matrix_index<size_t>;

//!\brief An offset type to move a matrix iterator in two-dimensional space.
//!\ingroup alignment_matrix
using matrix_offset = matrix_index<std::ptrdiff_t>;
} // namespace seqan3::detail
