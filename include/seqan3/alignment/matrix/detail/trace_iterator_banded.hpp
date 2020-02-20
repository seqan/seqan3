// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::trace_iterator_banded.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/matrix/detail/trace_iterator_base.hpp>

namespace seqan3::detail
{

/*!\brief A trace iterator for banded trace matrices.
 * \ingroup alignment_matrix
 *
 * \tparam matrix_iter_t The wrapped matrix iterator; must model seqan3::detail::two_dimensional_matrix_iterator and
 *                       the iterator's value type must be the same as seqan3::detail::trace_directions, i.e.
 *                       `std::same_as<std::iter_value_t<matrix_iter_t>, trace_directions>` must evaluate to `true`.
 *
 * \details
 *
 * This iterator follows a given trace in a banded trace matrix.
 */
template <two_dimensional_matrix_iterator matrix_iter_t>
class trace_iterator_banded : public trace_iterator_base<trace_iterator_banded<matrix_iter_t>, matrix_iter_t>
{
private:
    static_assert(std::same_as<std::iter_value_t<matrix_iter_t>, trace_directions>,
                  "Value type of the underlying iterator must be trace_directions.");

    //!\brief The type of the base class.
    using base_t = trace_iterator_base<trace_iterator_banded<matrix_iter_t>, matrix_iter_t>;

    //!\brief Befriend base class.
    friend base_t;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr trace_iterator_banded() = default; //!< Defaulted.
    constexpr trace_iterator_banded(trace_iterator_banded const &) = default; //!< Defaulted.
    constexpr trace_iterator_banded(trace_iterator_banded &&) = default; //!< Defaulted.
    constexpr trace_iterator_banded & operator=(trace_iterator_banded const &) = default; //!< Defaulted.
    constexpr trace_iterator_banded & operator=(trace_iterator_banded &&) = default; //!< Defaulted.
    ~trace_iterator_banded() = default; //!< Defaulted.

    /*!\brief Constructs from the underlying trace matrix iterator indicating the start of the trace path.
     * \param[in] matrix_iter The underlying matrix iterator.
     * \param[in] pivot_column The last column index which is still inside of the band in the first row of the
     *                         banded matrix.
     */
    template <typename index_t>
    constexpr trace_iterator_banded(matrix_iter_t const matrix_iter, column_index_type<index_t> const & pivot_column)
        noexcept :
        base_t{matrix_iter},
        pivot_column{static_cast<size_t>(pivot_column.get())}
    {}

    /*!\brief Constructs from the underlying trace matrix iterator indicating the start of the trace path.
     * \tparam other_matrix_iter_t The underlying matrix iterator type of `other`; the condition
     *                             `std::constructible_from<matrix_iter_t, other_matrix_iter_t>` must evaluate to `true`.
     * \param[in] other The underlying matrix iterator.
     *
     * \details
     *
     * Allows the conversion of non-const to const iterator.
     */
    template <two_dimensional_matrix_iterator other_matrix_iter_t>
    //!\cond
        requires std::constructible_from<matrix_iter_t, other_matrix_iter_t>
    //!\endcond
    constexpr trace_iterator_banded(trace_iterator_banded<other_matrix_iter_t> const other) noexcept : base_t{other}
    {}
    //!\}

    //!\copydoc seqan3::detail::trace_iterator_base::coordinate()
    [[nodiscard]] constexpr matrix_coordinate coordinate() const noexcept
    {
        auto coord = base_t::coordinate();
        coord.row += static_cast<int32_t>(coord.col - pivot_column);
        return coord;
    }

private:
    //!\copydoc seqan3::detail::trace_iterator_base::go_left
    constexpr void go_left(matrix_iter_t & iter) const noexcept
    {
        // Note, in the banded matrix, the columns are virtually shifted by one cell.
        // So going left means go to the previous column and then one row down.
        iter -= matrix_offset{row_index_type{-1}, column_index_type{1}};
    }

    //!\copydoc seqan3::detail::trace_iterator_base::go_up
    constexpr void go_diagonal(matrix_iter_t & iter) const noexcept
    {
        // Note, in the banded matrix, the columns are virtually shifted by one cell.
        // So going diagonal means go to the previous column and stay in the same row.
        iter -= matrix_offset{row_index_type{0}, column_index_type{1}};
    }

    size_t pivot_column{}; //!< The largest column index which is inside of the band in the first row of the matrix.
};

} // namespace seqan3::detail
