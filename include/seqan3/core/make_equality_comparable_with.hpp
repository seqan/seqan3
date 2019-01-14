// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::make_equality_comparable_with.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{

/*!\brief Adds comparison operators to the derived type.
 * \ingroup core
 * \tparam derived_t    The type of the class that inherits from this class.
 * \tparam comparison_t The type that is made comparable with `derived_t`.
 *
 * \details
 *
 * This CRTP-base class adds comparison operators to `derived_t`, i.e. `derived_t` can be compared with itself
 * and with the specified second template parameter. This base class simplifies the implementation of all necessary
 * comparison operators if the comparison type cannot be converted into the type `derived_t` due to missing
 * conversion constructors in `derived_t` or conversion operators in the comparison type. A typical example for
 * this could be proxy types returned by an iterator, which represents another value, but cannot itself be constructed
 * from it.
 */
template <typename derived_t, typename comparison_t>
//!\cond
    requires std::EqualityComparable<comparison_t>
//!\endcond
class make_equality_comparable_with
{
private:

    //!\brief Befriend the derived class.
    friend derived_t;

    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    make_equality_comparable_with() = default;
    make_equality_comparable_with(make_equality_comparable_with const &) = default;
    make_equality_comparable_with(make_equality_comparable_with &&) = default;
    make_equality_comparable_with & operator=(make_equality_comparable_with const &) = default;
    make_equality_comparable_with & operator=(make_equality_comparable_with &&) = default;
    ~make_equality_comparable_with() = default;
    //!}

    //!\brief Casts `this` to the derived type.
    constexpr derived_t & underlying()
    {
        return static_cast<derived_t &>(*this);
    }

    //!\brief Casts `this` to the derived type.
    constexpr derived_t const & underlying() const
    {
        return static_cast<derived_t const &>(*this);
    }

public:

    /*!\name Comparison operators
     * \{
     */

    /*!\brief Tests for equality.
     * \param[in] rhs The right operand of the comparison.
     */
    constexpr bool operator==(comparison_t const & rhs) const noexcept
    {
        return underlying().compare_value() == rhs;
    }

    //!\copydoc operator==
    constexpr bool operator==(make_equality_comparable_with const & rhs) const noexcept
    {
        return *this == rhs.underlying().compare_value();
    }

    /*!\brief Tests for equality.
     * \param[in] lhs The left operand of the comparison.
     * \param[in] rhs The right operand of the comparison.
     *
     * \details
     *
     * Free function to allow comparison with `make_equality_comparable_with` on the right-hand side and
     * `comparison_t` on the left-hand side.
     */
    friend constexpr bool operator==(comparison_t const & lhs,
                                     make_equality_comparable_with const & rhs) noexcept
    {
        return rhs == lhs;
    }

    /*!\brief Tests for inequality.
     * \param[in] rhs The right operand of the comparison.
     */
    constexpr bool operator!=(comparison_t const & rhs) const noexcept
    {
        return underlying().compare_value() != rhs;
    }

    //!\copydoc operator!=
    constexpr bool operator!=(make_equality_comparable_with const & rhs) const noexcept
    {
        return *this != rhs.underlying().compare_value();
    }

    /*!\brief Tests for inequality.
     * \param[in] lhs The left operand of the comparison.
     * \param[in] rhs The right operand of the comparison.
     *
     * \details
     *
     * Free function to allow comparison with `make_equality_comparable_with` on the right-hand side and
     * `comparison_t` on the left-hand side.
     */
    friend constexpr bool operator!=(comparison_t const & lhs,
                                     make_equality_comparable_with const & rhs) noexcept
    {
        return rhs != lhs;
    }

    /*!\brief Tests for less than.
     * \param[in] rhs The right operand of the comparison.
     */
    constexpr bool operator<(comparison_t const & rhs) const noexcept
    {
        return underlying().compare_value() < rhs;
    }

    //!\copydoc operator<
    constexpr bool operator<(make_equality_comparable_with const & rhs) const noexcept
    {
        return *this < rhs.underlying().compare_value();
    }

    /*!\brief Tests for less than.
     * \param[in] lhs The left operand of the comparison.
     * \param[in] rhs The right operand of the comparison.
     *
     * \details
     *
     * Free function to allow comparison with `make_equality_comparable_with` on the right-hand side and
     * `comparison_t` on the left-hand side.
     */
    friend constexpr bool operator<(comparison_t const & lhs,
                                    make_equality_comparable_with const & rhs) noexcept
    {
        return rhs >= lhs;
    }

    /*!\brief Tests for less than or equality.
     * \param[in] rhs The right operand of the comparison.
     */
    constexpr bool operator<=(comparison_t const & rhs) const noexcept
    {
        return underlying().compare_value() <= rhs;
    }

    //!\copydoc operator<=
    constexpr bool operator<=(make_equality_comparable_with const & rhs) const noexcept
    {
        return *this <= rhs.underlying().compare_value();
    }

    /*!\brief Tests for less than or equality.
     * \param[in] lhs The left operand of the comparison.
     * \param[in] rhs The right operand of the comparison.
     *
     * \details
     *
     * Free function to allow comparison with `make_equality_comparable_with` on the right-hand side and
     * `comparison_t` on the left-hand side.
     */
    friend constexpr bool operator<=(comparison_t const & lhs,
                                     make_equality_comparable_with const & rhs) noexcept
    {
        return rhs >= lhs;
    }

    /*!\brief Tests for greater than.
     * \param[in] rhs The right operand of the comparison.
     */
    constexpr bool operator>(comparison_t const & rhs) const noexcept
    {
        return underlying().compare_value() > rhs;
    }

    //!\copydoc operator>
    constexpr bool operator>(make_equality_comparable_with const & rhs) const noexcept
    {
        return *this > rhs.underlying().compare_value();
    }

    /*!\brief Tests for greater than.
     * \param[in] lhs The left operand of the comparison.
     * \param[in] rhs The right operand of the comparison.
     *
     * \details
     *
     * Free function to allow comparison with `make_equality_comparable_with` on the right-hand side and
     * `comparison_t` on the left-hand side.
     */
    friend constexpr bool operator>(comparison_t const & lhs,
                                    make_equality_comparable_with const & rhs) noexcept
    {
        return rhs <= lhs;
    }

    /*!\brief Tests for greater than or equality.
     * \param[in] rhs The right operand of the comparison.
     */
    constexpr bool operator>=(comparison_t const & rhs) const noexcept
    {
        return underlying().compare_value() >= rhs;
    }

    //!\copydoc operator>=
    constexpr bool operator>=(make_equality_comparable_with const & rhs) const noexcept
    {
        return *this >= rhs.underlying().compare_value();
    }

    /*!\brief Tests for greater than or equality.
     * \param[in] lhs The left operand of the comparison.
     * \param[in] rhs The right operand of the comparison.
     *
     * \details
     *
     * Free function to allow comparison with `make_equality_comparable_with` on the right-hand side and
     * `comparison_t` on the left-hand side.
     */
    friend constexpr bool operator>=(comparison_t const & lhs,
                                     make_equality_comparable_with const & rhs) noexcept
    {
        return rhs <= lhs;
    }
    //!\}
};
} // namespace seqan3
