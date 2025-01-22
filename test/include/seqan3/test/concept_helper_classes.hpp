// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides helper classes for testing concepts.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

// Helper classes for testing concepts.
class type_b;
class type_c;
class type_d;

/*!\brief Helper class for testing concepts.
 * \details
 * * weakly_ordered with itself
 * * totally_ordered with type_b
 */
class type_a
{
public:
    /*!\name Default constructors and assignments.
     * \{
     */
    type_a() = default;                           //!< Defaulted.
    type_a(type_a const &) = default;             //!< Defaulted.
    type_a & operator=(type_a const &) = default; //!< Defaulted.
    type_a(type_a &&) = default;                  //!< Defaulted.
    type_a & operator=(type_a &&) = default;      //!< Defaulted.
    ~type_a() = default;                          //!< Defaulted.
    //!\}

    /*!\name Comparison operators.
     * \{
     */
    //!@{ Compares to another type.
    friend constexpr bool operator<(type_a const & lhs, type_a const & rhs);
    friend constexpr bool operator>(type_a const & lhs, type_a const & rhs);
    friend constexpr bool operator<=(type_a const & lhs, type_a const & rhs);
    friend constexpr bool operator>=(type_a const & lhs, type_a const & rhs);

    friend constexpr bool operator==(type_a const & lhs, type_b const & rhs);
    friend constexpr bool operator!=(type_a const & lhs, type_b const & rhs);
    friend constexpr bool operator<(type_a const & lhs, type_b const & rhs);
    friend constexpr bool operator>(type_a const & lhs, type_b const & rhs);
    friend constexpr bool operator<=(type_a const & lhs, type_b const & rhs);
    friend constexpr bool operator>=(type_a const & lhs, type_b const & rhs);
    //!@}
    //!\}
};

/*!\brief Move-only helper class for testing concepts.
 * \details
 * * totally_ordered with itself, type_a, type_d
 */
class type_b : public type_a
{
public:
    /*!\name Default constructors and assignments.
     * \{
     */
    type_b() = default;                          //!< Defaulted.
    type_b(type_b const &) = delete;             //!< Deleted.
    type_b & operator=(type_b const &) = delete; //!< Deleted.
    type_b(type_b &&) = default;                 //!< Defaulted.
    type_b & operator=(type_b &&) = default;     //!< Defaulted.
    ~type_b() = default;                         //!< Defaulted.
    //!\}

    /*!\name Comparison operators.
     * \{
     */
    //!@{ Compares to another type.
    friend constexpr bool operator==(type_b const & lhs, type_b const & rhs);
    friend constexpr bool operator!=(type_b const & lhs, type_b const & rhs);
    friend constexpr bool operator<(type_b const & lhs, type_b const & rhs);
    friend constexpr bool operator>(type_b const & lhs, type_b const & rhs);
    friend constexpr bool operator<=(type_b const & lhs, type_b const & rhs);
    friend constexpr bool operator>=(type_b const & lhs, type_b const & rhs);

    friend constexpr bool operator==(type_b const & lhs, type_d const & rhs);
    friend constexpr bool operator!=(type_b const & lhs, type_d const & rhs);
    friend constexpr bool operator<(type_b const & lhs, type_d const & rhs);
    friend constexpr bool operator>(type_b const & lhs, type_d const & rhs);
    friend constexpr bool operator<=(type_b const & lhs, type_d const & rhs);
    friend constexpr bool operator>=(type_b const & lhs, type_d const & rhs);
    //!@}
    //!\}
};

/*!\brief Helper class with custom constructors for testing concepts.
 * \details
 * * equality_comparable with itself
 */
class type_c
{
public:
    /*!\name Default constructors and assignments.
     * \{
     */
    type_c() = default;                           //!< Defaulted.
    type_c(type_c const &) = default;             //!< Defaulted.
    type_c & operator=(type_c const &) = default; //!< Defaulted.
    type_c(type_c &&) = default;                  //!< Defaulted.
    type_c & operator=(type_c &&) = default;      //!< Defaulted.
    ~type_c() = default;                          //!< Defaulted.
    //!\}

    //!\brief Construct from type_b.
    type_c(type_b const &)
    {}

    //!\brief Construct from type_a.
    explicit type_c(type_a const &)
    {}

    /*!\name Comparison operators.
     * \{
     */
    //!@{ Compares to another type.
    friend constexpr bool operator==(type_c const & lhs, type_c const & rhs);
    friend constexpr bool operator!=(type_c const & lhs, type_c const & rhs);
    //!@}
    //!\}
};

/*!\brief Unconstructible helper class for testing concepts.
 * \details
 * * totally_ordered with itself, type_b
 */
class type_d : public type_b
{
public:
    /*!\name Default constructors and assignments.
     * \{
     */
    type_d() = delete;                           //!< Deleted.
    type_d(type_d const &) = delete;             //!< Deleted.
    type_d & operator=(type_d const &) = delete; //!< Deleted.
    type_d(type_d &&) = delete;                  //!< Deleted.
    type_d & operator=(type_d &&) = delete;      //!< Deleted.
    ~type_d() = delete;                          //!< Deleted.
    //!\}

    /*!\name Comparison operators.
     * \{
     */
    //!@{ Compares to another type.
    friend constexpr bool operator==(type_d const & lhs, type_d const & rhs);
    friend constexpr bool operator!=(type_d const & lhs, type_d const & rhs);
    friend constexpr bool operator<(type_d const & lhs, type_d const & rhs);
    friend constexpr bool operator>(type_d const & lhs, type_d const & rhs);
    friend constexpr bool operator<=(type_d const & lhs, type_d const & rhs);
    friend constexpr bool operator>=(type_d const & lhs, type_d const & rhs);
    //!@}
    //!\}
};
