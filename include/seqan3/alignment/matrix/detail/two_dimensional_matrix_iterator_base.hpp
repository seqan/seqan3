// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::two_dimensional_matrix_iterator_base.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <ranges>
#include <type_traits>

#include <seqan3/alignment/matrix/detail/matrix_coordinate.hpp>
#include <seqan3/core/detail/iterator_traits.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

namespace seqan3::detail
{

/*!\brief Selects the major order of the matrix.
 * \ingroup alignment_matrix
 *
 * \details
 *
 * This enum is used to select between column and row major order access patterns for
 * seqan3::detail::two_dimensional_matrix.
 * This matrix type stores a two-dimensional matrix in a flattened one-dimensional vector, whose access orientation
 * can be adapted using this policy.
 *
 * \see seqan3::detail::two_dimensional_matrix_iterator_base
 */
enum struct matrix_major_order : uint8_t
{
    column, //!< Accesses matrix in column major order.
    row     //!< Accesses matrix in row major order.
};

/*!\brief A crtp-base class for iterators over seqan3::detail::two_dimensional_matrix.
 * \ingroup alignment_matrix
 * \implements seqan3::detail::two_dimensional_matrix_iterator
 * \tparam derived_t The derived type that implements this base iterator.
 * \tparam order The seqan3::detail::matrix_major_order to use.
 *
 * \details
 *
 * \note This is a crtp_base class and cannot be instantiated directly. You need to define an iterator type that
 *       implements this abstract base class.
 *
 * This iterator provides a two-dimensional access interface over the seqan3::detail::two_dimensional_matrix, which
 * stores the values in a one-dimensional vector.
 * In addition to the regular interface using seqan3::detail::two_dimensional_matrix_iterator_base::difference_type<dummy_t>
 * this iterator offers special operators for advancing the iterator in a two-dimensional layout using
 * seqan3::detail::matrix_offset. The underlying iterator is moved along the respective row and column offset
 * according to the specialised seqan3::detail::matrix_major_order.
 *
 * The regular interface moves the wrapped iterator according to the specified seqan3::detail::matrix_major_order, i.e.
 * if `order` is seqan3::detail::matrix_major_order::column it advances the iterator first in vertical dimension and
 * second in horizontal dimension and if `order` is seqan3::detail::matrix_major_order::row the other way around.
 *
 * ### Requirements of the derived type
 *
 * Types that implement this abstract base class must provide the associated types for iterators and
 * shall implement the following functions in order to work properly:
 *
 * * seqan3::detail::two_dimensional_matrix_iterator_base::operator+=(matrix_offset)
 * * seqan3::detail::two_dimensional_matrix_iterator_base::coordinate
 *
 */
template <typename derived_t, matrix_major_order order>
class two_dimensional_matrix_iterator_base
{
private:
    //!\brief Befriend the derived type.
    friend derived_t;

    //!\brief Befriend other base class types for const iterator compatibility.
    template <typename other_derived_t, matrix_major_order other_order>
    friend class two_dimensional_matrix_iterator_base;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Defaulted.
    constexpr two_dimensional_matrix_iterator_base() = default;
    //!\brief Defaulted.
    constexpr two_dimensional_matrix_iterator_base(two_dimensional_matrix_iterator_base const &) = default;
    //!\brief Defaulted.
    constexpr two_dimensional_matrix_iterator_base(two_dimensional_matrix_iterator_base &&) = default;
    //!\brief Defaulted.
    constexpr two_dimensional_matrix_iterator_base & operator=(two_dimensional_matrix_iterator_base const &) = default;
    //!\brief Defaulted.
    constexpr two_dimensional_matrix_iterator_base & operator=(two_dimensional_matrix_iterator_base &&) = default;
    //!\brief Defaulted.
    ~two_dimensional_matrix_iterator_base() = default;
    //!\}

    //!\brief Helper template definition to get the difference type of the derived type.
    template <typename _derived_t>
    using difference_type = typename _derived_t::difference_type;

    //!\brief Helper template definition to get the reference type of the derived type.
    template <typename _derived_t>
    using reference = typename _derived_t::reference;

    //!\brief Helper template definition to get the pointer type of the derived type.
    template <typename _derived_t>
    using pointer = typename _derived_t::pointer;

public:
    /*!\name Element access
     * \{
     */
    //!\brief Returns a reference to the pointed to element.
    template <typename dummy_t = derived_t>
    constexpr reference<dummy_t> operator*() const noexcept
    {
        return *as_derived().host_iter;
    }

    //!\brief Returns a reference to the pointed-to-element after advancing the iterator by the given offset.
    template <typename dummy_t = derived_t>
    constexpr reference<dummy_t> operator[](std::iter_difference_t<dummy_t> const offset) const noexcept
    {
        return *(as_derived() + offset);
    }

    //!\brief Returns a reference to the pointed-to-element after advancing the iterator by the given offset.
    template <typename dummy_t = derived_t>
    constexpr reference<dummy_t> operator[](matrix_offset const & offset) const noexcept
    {
        return *(as_derived() + offset);
    }

    //!\brief Returns a pointer to the pointed-to-element.
    template <typename dummy_t = derived_t>
    constexpr pointer<dummy_t> operator->() const noexcept
    {
        return std::addressof(*as_derived().host_iter);
    }

    /*!\brief Returns the current position of the iterator as a seqan3::detail::matrix_coordinate
     *
     * \details
     *
     * The position of the iterator is stored as a seqan3::detail::matrix_coordinate mapping the one-dimensional
     * vector position to a two-dimensional point coordinate.
     *
     * \note This function shall be implented in the derived type.
     */

    SEQAN3_DOXYGEN_ONLY(constexpr seqan3::detail::matrix_coordinate coordinate() const noexcept {})
    //!\}

    /*!\name Arithmetic operators
     * \{
     */

    /*!\brief Advances the iterator by the given offset in the respective matrix dimensions.
     * \param[in] offset The matrix offset used to advance the iterator.
     * \returns The advanced iterator.
     *
     * \details
     *
     * Advances the given host iterator (the iterator over the one-dimensional vector) by the given matrix coordinate
     * offset. Independent of the seqan3::detail::matrix_major_order the
     * host iterator is advanced to the correct position within the underlying one-dimensional vector as if it where
     * a two-dimensional matrix.
     *
     * \note This operator shall be implemented by the derived type and all other arithmetic operators will delegate
     *       to this operator.
     */
    SEQAN3_DOXYGEN_ONLY(constexpr derived_t & operator+=(matrix_offset const & offset) noexcept {})

    //!\brief Advances the iterator by one following the given matrix major order.
    constexpr derived_t & operator++() noexcept
    {
        if constexpr (order == matrix_major_order::column)
            return as_derived() += matrix_offset{row_index_type{1}, column_index_type{0}};
        else
            return as_derived() += matrix_offset{row_index_type{0}, column_index_type{1}};
    }

    //!\brief Returns an iterator incremented by one following the given matrix major order.
    constexpr derived_t operator++(int) noexcept
    {
        derived_t previous{as_derived()};
        ++(*this);
        return previous;
    }

    //!\brief Advances the iterator by `offset` following the given matrix major order.
    template <typename dummy_t = derived_t>
    constexpr derived_t & operator+=(std::iter_difference_t<dummy_t> const offset) noexcept
    {
        if constexpr (order == matrix_major_order::column)
            return as_derived() += matrix_offset{row_index_type{offset}, column_index_type{0}};
        else
            return as_derived() += matrix_offset{row_index_type{0}, column_index_type{offset}};
    }

    //!\brief Returns an iterator advanced by `offset` following the given matrix major order.
    template <typename dummy_t = derived_t>
    constexpr derived_t operator+(std::iter_difference_t<dummy_t> const offset) const noexcept
    {
        derived_t next{as_derived()};
        next += offset;
        return next;
    }

    //!\brief Returns an iterator advanced by `offset` following the given matrix major order.
    template <typename dummy_t = derived_t>
    constexpr friend derived_t operator+(std::iter_difference_t<dummy_t> const offset, derived_t const iter)
    {
        return iter + offset;
    }

    //!\brief Returns an iterator advanced by `offset` in the respective dimensions.
    constexpr derived_t operator+(matrix_offset const & offset) const noexcept
    {
        derived_t next{as_derived()};
        next += offset;
        return next;
    }

    //!\brief Returns an iterator advanced by `offset` in the respective dimensions.
    constexpr friend derived_t operator+(matrix_offset const & offset, derived_t const iter)
    {
        return iter + offset;
    }

    //!\brief Advances the iterator by minus one following the given matrix major order.
    constexpr derived_t & operator--() noexcept
    {
        if constexpr (order == matrix_major_order::column)
            return as_derived() += matrix_offset{row_index_type{-1}, column_index_type{0}};
        else
            return as_derived() += matrix_offset{row_index_type{0}, column_index_type{-1}};
    }

    //!\brief Returns an iterator decremented by one following the given matrix major order.
    constexpr derived_t operator--(int) noexcept
    {
        derived_t previous{as_derived()};
        --(*this);
        return previous;
    }

    //!\brief Advances the iterator by `offset` following the given matrix major order.
    template <typename dummy_t = derived_t>
    constexpr derived_t & operator-=(std::iter_difference_t<dummy_t> const offset) noexcept
    {
        return *this += -offset;
    }

    //!\brief Returns an iterator advanced by `offset` following the given matrix major order.
    template <typename dummy_t = derived_t>
    constexpr derived_t operator-(std::iter_difference_t<dummy_t> const offset) const noexcept
    {
        derived_t next{as_derived()};
        next -= offset;
        return next;
    }

    //!\brief Returns an iterator advanced by `offset` in the respective dimensions.
    constexpr derived_t & operator-=(matrix_offset const & offset) noexcept
    {
        return as_derived() += matrix_offset{row_index_type{-offset.row}, column_index_type{-offset.col}};
    }

    //!\brief Returns an iterator advanced by `offset` in the respective dimensions.
    constexpr derived_t operator-(matrix_offset const & offset) const noexcept
    {
        derived_t next{as_derived()};
        next -= offset;
        return next;
    }

    //!\brief Returns the distance between two iterators.
    template <typename dummy_t = derived_t>
    friend constexpr std::iter_difference_t<dummy_t> operator-(derived_t const lhs, derived_t const rhs) noexcept
    {
        return lhs.as_host_iter() - rhs.as_host_iter();
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Returns `true` if both iterators are equal, `false` otherwise.
    // What if the derived type is different?
    // How can make sure this is the same type?
    template <typename other_derived_t>
        requires std::constructible_from<derived_t, other_derived_t>
              || std::constructible_from<other_derived_t, derived_t>
    constexpr bool operator==(two_dimensional_matrix_iterator_base<other_derived_t, order> const & rhs) const noexcept
    {
        return as_derived().host_iter == rhs.as_derived().host_iter;
    }

    //!\brief Returns `true` if both iterators are unequal, `false` otherwise.
    template <typename other_derived_t>
        requires std::constructible_from<derived_t, other_derived_t>
              || std::constructible_from<other_derived_t, derived_t>
    constexpr bool operator!=(two_dimensional_matrix_iterator_base<other_derived_t, order> const & rhs) const noexcept
    {
        return !(*this == rhs);
    }

    //!\brief Checks if `lhs` is smaller than `rhs`.
    template <typename other_derived_t>
        requires std::constructible_from<derived_t, other_derived_t>
              || std::constructible_from<other_derived_t, derived_t>
    constexpr bool operator<(two_dimensional_matrix_iterator_base<other_derived_t, order> const & rhs) const noexcept
    {
        return as_derived().host_iter < rhs.as_derived().host_iter;
    }

    //!\brief Checks if `lhs` is smaller than or equal to `rhs`.
    template <typename other_derived_t>
        requires std::constructible_from<derived_t, other_derived_t>
              || std::constructible_from<other_derived_t, derived_t>
    constexpr bool operator<=(two_dimensional_matrix_iterator_base<other_derived_t, order> const & rhs) const noexcept
    {
        return as_derived().host_iter <= rhs.as_derived().host_iter;
    }

    //!\brief Checks if `lhs` is greater than `rhs`.
    template <typename other_derived_t>
        requires std::constructible_from<derived_t, other_derived_t>
              || std::constructible_from<other_derived_t, derived_t>
    constexpr bool operator>(two_dimensional_matrix_iterator_base<other_derived_t, order> const & rhs) const noexcept
    {
        return as_derived().host_iter > rhs.as_derived().host_iter;
    }

    //!\brief Checks if `lhs` is greater than or equal to `rhs`.
    template <typename other_derived_t>
        requires std::constructible_from<derived_t, other_derived_t>
              || std::constructible_from<other_derived_t, derived_t>
    constexpr bool operator>=(two_dimensional_matrix_iterator_base<other_derived_t, order> const & rhs) const noexcept
    {
        return as_derived().host_iter >= rhs.as_derived().host_iter;
    }
    //!\}

private:
    //!\brief Return the host_iter of the derived type.
    constexpr auto const & as_host_iter() const
    {
        return as_derived().host_iter;
    }

    //!\brief Cast this to derived type.
    constexpr derived_t & as_derived()
    {
        return static_cast<derived_t &>(*this);
    }

    //!\copydoc seqan3::detail::two_dimensional_matrix_iterator_base::as_derived
    constexpr derived_t const & as_derived() const
    {
        return static_cast<derived_t const &>(*this);
    }

    // matrix_iterator_t host_iter{};    //!< The wrapped matrix iterator.
};
} // namespace seqan3::detail
