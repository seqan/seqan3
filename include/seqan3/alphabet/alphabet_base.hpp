// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::alphabet_base.
 */

#pragma once

#include <cassert>
#include <concepts>
#include <type_traits>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/utility/detail/integer_traits.hpp>

namespace seqan3
{

/*!\brief A CRTP-base that makes defining a custom alphabet easier.
 * \ingroup alphabet
 * \tparam derived_type The CRTP parameter type.
 * \tparam size         The size of the alphabet.
 * \tparam char_t       The character type of the alphabet (set this to `void` when defining just a
 *                      seqan3::semialphabet).
 *
 * \details
 *
 * You can use this class to define your own alphabet, but types are not required to be based on it to model
 * seqan3::alphabet, it is purely a way to avoid code duplication.
 *
 * The base class represents the alphabet value as the rank and automatically deduces the rank type from the size and
 * defines all required member functions. The derived type needs to define only the following two static
 * member functions (can be private if the base class is befriended):
 *
 *   * `static constexpr char_type rank_to_char(rank_type const rank);` that defines for every possible rank value the
 *     corresponding char value.
 *     (The implementation can be a lookup-table or an arithmetic expression.)
 *   * `static constexpr rank_type char_to_rank(char_type const chr);` that defines for every possible character value
 *     the corresponding rank value (adapt size if char_type isn't `char`).
 *     (The implementation can be a lookup-table or an arithmetic expression.)
 *
 * ### Example
 *
 * This creates an alphabet called `ab` which has size two and the two letters 'A' and 'B':
 * \include test/snippet/alphabet/alphabet_base.cpp
 *
 * \stableapi{Since version 3.1.}
 */
template <typename derived_type, size_t size, typename char_t = char>
class alphabet_base
{
protected:
    static_assert(size != 0, "alphabet size must be >= 1"); // == 1 is handled below in separate specialisation

    /*!\name Member types
     * \{
     */
    /*!\brief The char representation; conditional needed to make semi alphabet definitions legal.
     * \details We need a return type for seqan3::alphabet_base::to_char and seqan3::alphabet_base::assign_char other
     *          than void to make these in-class definitions valid when `char_t` is void.
     *
     * \attention Please use seqan3::alphabet_char_t to access this type.
     *
     * \stableapi{Since version 3.1.}
     */
    using char_type = std::conditional_t<std::same_as<char_t, void>, char, char_t>;

    /*!\brief The type of the alphabet when represented as a number (e.g. via to_rank()).
     *
     * \attention Please use seqan3::alphabet_rank_t to access this type.
     *
     * \stableapi{Since version 3.1.}
     */
    using rank_type = detail::min_viable_uint_t<size - 1>;
    //!\}

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr alphabet_base() noexcept = default;                                  //!< Defaulted.
    constexpr alphabet_base(alphabet_base const &) noexcept = default;             //!< Defaulted.
    constexpr alphabet_base(alphabet_base &&) noexcept = default;                  //!< Defaulted.
    constexpr alphabet_base & operator=(alphabet_base const &) noexcept = default; //!< Defaulted.
    constexpr alphabet_base & operator=(alphabet_base &&) noexcept = default;      //!< Defaulted.
    ~alphabet_base() noexcept = default;                                           //!< Defaulted.

    //!\}

    /*!\name Read functions
     * \{
     */
    /*!\brief Return the letter as a character of char_type.
     *
     * \details
     *
     * Provides an implementation for seqan3::to_char, required to model seqan3::alphabet.
     *
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
     *
     * Guaranteed not to throw.
     *
     * \stableapi{Since version 3.1.}
     */
    constexpr char_type to_char() const noexcept
        requires (!std::same_as<char_t, void>)
    {
        return derived_type::rank_to_char(rank);
    }

    /*!\brief Return the letter's numeric value (rank in the alphabet).
     *
     * \details
     *
     * Provides an implementation for seqan3::to_rank, required to model seqan3::semialphabet.
     *
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
     *
     * Guaranteed not to throw.
     *
     * \stableapi{Since version 3.1.}
     */
    constexpr rank_type to_rank() const noexcept
    {
        return rank;
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    /*!\brief Assign from a character, implicitly converts invalid characters.
     * \param chr The character to be assigned.
     *
     * \details
     *
     * Provides an implementation for seqan3::assign_char_to, required to model seqan3::alphabet.
     *
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
     *
     * Guaranteed not to throw.
     *
     * \stableapi{Since version 3.1.}
     */
    constexpr derived_type & assign_char(char_type const chr) noexcept
        requires (!std::same_as<char_t, void>)
    {
        rank = derived_type::char_to_rank(chr);
        return static_cast<derived_type &>(*this);
    }

    /*!\brief Assign from a numeric value.
     * \param c The rank to be assigned.
     *
     * \details
     *
     * Provides an implementation for seqan3::assign_rank_to, required to model seqan3::semialphabet.
     *
     * ###Complexity
     *
     * Constant.
     *
     * ###Exceptions
     *
     * Guaranteed not to throw.
     *
     * \stableapi{Since version 3.1.}
     */
    constexpr derived_type & assign_rank(rank_type const c) noexcept
    {
        assert(static_cast<size_t>(c) < static_cast<size_t>(alphabet_size));
        rank = c;
        return static_cast<derived_type &>(*this);
    }
    //!\}

    /*!\brief The size of the alphabet, i.e. the number of different values it can take.
     *
     * \stableapi{Since version 3.1.}
     */
    static constexpr detail::min_viable_uint_t<size> alphabet_size = size;

    //!\name Comparison operators
    //!\{

    /*!\brief Checks whether the letters `lhs` and `rhs` are equal.
     *
     * \stableapi{Since version 3.1.}
     */
    friend constexpr bool operator==(derived_type const lhs, derived_type const rhs) noexcept
    {
        return seqan3::to_rank(lhs) == seqan3::to_rank(rhs);
    }

    /*!\brief Checks whether the letters `lhs` and `rhs` are unequal.
     *
     * \stableapi{Since version 3.1.}
     */
    friend constexpr bool operator!=(derived_type const lhs, derived_type const rhs) noexcept
    {
        return seqan3::to_rank(lhs) != seqan3::to_rank(rhs);
    }

    /*!\brief Checks whether the letter `lhs` is smaller than `rhs`.
     *
     * \stableapi{Since version 3.1.}
     */
    friend constexpr bool operator<(derived_type const lhs, derived_type const rhs) noexcept
    {
        return seqan3::to_rank(lhs) < seqan3::to_rank(rhs);
    }

    /*!\brief Checks whether the letter `lhs` is greater than `rhs`.
     *
     * \stableapi{Since version 3.1.}
     */
    friend constexpr bool operator>(derived_type const lhs, derived_type const rhs) noexcept
    {
        return seqan3::to_rank(lhs) > seqan3::to_rank(rhs);
    }

    /*!\brief Checks whether the letter `lhs` is smaller than or equal to `rhs`.
     *
     * \stableapi{Since version 3.1.}
     */
    friend constexpr bool operator<=(derived_type const lhs, derived_type const rhs) noexcept
    {
        return seqan3::to_rank(lhs) <= seqan3::to_rank(rhs);
    }

    /*!\brief Checks whether the letter `lhs` is bigger than or equal to `rhs`.
     *
     * \stableapi{Since version 3.1.}
     */
    friend constexpr bool operator>=(derived_type const lhs, derived_type const rhs) noexcept
    {
        return seqan3::to_rank(lhs) >= seqan3::to_rank(rhs);
    }
    //!\}

private:
    //!\brief The value of the alphabet letter is stored as the rank.
    rank_type rank;
};

} // namespace seqan3
