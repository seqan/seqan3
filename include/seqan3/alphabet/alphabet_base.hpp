// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::alphabet_base.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/int_types.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/type_traits>

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
 * defines all required member functions. The derived type needs to define only the following two tables as static
 * member variables (can be private if the base class is befriended):
 *
 *   * `static std::array<char_type, alphabet_size> constexpr rank_to_char` that defines for every possible rank value
 *     the corresponding char value.
 *   * `static std::array<rank_type, 256> constexpr char_to_rank` that defines for every possible character value the
 *     corresponding rank value (adapt size if char_type isn't `char`).
 *
 * ### Example
 *
 * This creates an alphabet called `ab` which has size two and the two letters 'A' and 'B':
 * \include test/snippet/alphabet/detail/alphabet_base.cpp
 *
 */
template <typename derived_type, size_t size, typename char_t = char>
class alphabet_base
{
protected:
    static_assert(size != 0, "alphabet size must be >= 1"); // == 1 is handled below in separate specialisation

    /*!\name Member types
     * \{
     */
    //!\brief The char representation; conditional needed to make semi alphabet definitions legal.
    using char_type = std::conditional_t<std::same_as<char_t, void>, char, char_t>;
    //!\brief The type of the alphabet when represented as a number (e.g. via to_rank()).
    using rank_type = detail::min_viable_uint_t<size - 1>;
    //!\}

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr alphabet_base()                                   noexcept = default; //!< Defaulted.
    constexpr alphabet_base(alphabet_base const &)              noexcept = default; //!< Defaulted.
    constexpr alphabet_base(alphabet_base &&)                   noexcept = default; //!< Defaulted.
    constexpr alphabet_base & operator=(alphabet_base const &)  noexcept = default; //!< Defaulted.
    constexpr alphabet_base & operator=(alphabet_base &&)       noexcept = default; //!< Defaulted.
    ~alphabet_base()                                            noexcept = default; //!< Defaulted.
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
     */
    constexpr char_type to_char() const noexcept
    //!\cond
        requires (!std::same_as<char_t, void>)
    //!\endcond
    {
        return derived_type::rank_to_char[rank];
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
     * \param c The character to be assigned.
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
     */
    constexpr derived_type & assign_char(char_type const c) noexcept
    //!\cond
        requires (!std::same_as<char_t, void>)
    //!\endcond
    {
        using index_t = std::make_unsigned_t<char_type>;
        rank = derived_type::char_to_rank[static_cast<index_t>(c)];
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
     */
    constexpr derived_type & assign_rank(rank_type const c) noexcept
    {
        assert(static_cast<size_t>(c) < static_cast<size_t>(alphabet_size));
        rank = c;
        return static_cast<derived_type &>(*this);
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static detail::min_viable_uint_t<size> constexpr alphabet_size = size;

    //!\name Comparison operators
    //!\{

    //!\brief Checks whether the letters `lhs` and `rhs` are equal.
    friend constexpr bool operator==(derived_type const lhs, derived_type const rhs) noexcept
    {
        return seqan3::to_rank(lhs) == seqan3::to_rank(rhs);
    }

    //!\brief Checks whether the letters `lhs` and `rhs` are unequal.
    friend constexpr bool operator!=(derived_type const lhs, derived_type const rhs) noexcept
    {
        return seqan3::to_rank(lhs) != seqan3::to_rank(rhs);
    }

    //!\brief Checks whether the letter `lhs` is smaller than `rhs`.
    friend constexpr bool operator<(derived_type const lhs, derived_type const rhs) noexcept
    {
        return seqan3::to_rank(lhs) < seqan3::to_rank(rhs);
    }

    //!\brief Checks whether the letter `lhs` is greater than `rhs`.
    friend constexpr bool operator>(derived_type const lhs, derived_type const rhs) noexcept
    {
        return seqan3::to_rank(lhs) > seqan3::to_rank(rhs);
    }

    //!\brief Checks whether the letter `lhs` is smaller than or equal to `rhs`.
    friend constexpr bool operator<=(derived_type const lhs, derived_type const rhs) noexcept
    {
        return seqan3::to_rank(lhs) <= seqan3::to_rank(rhs);
    }

    //!\brief Checks whether the letter `lhs` is bigger than or equal to `rhs`.
    friend constexpr bool operator>=(derived_type const lhs, derived_type const rhs) noexcept
    {
        return seqan3::to_rank(lhs) >= seqan3::to_rank(rhs);
    }
    //!\}

private:
    //!\brief The value of the alphabet letter is stored as the rank.
    rank_type rank{};
};

/*!\brief Specialisation of seqan3::alphabet_base for alphabets of size 1.
 * \ingroup alphabet
 * \tparam derived_type The CRTP parameter type.
 * \tparam char_t The character type (always set to `char` for alphabets of size 1 and to `void` for semi alphabets of
 *                size 1).
 *
 * \details
 *
 * This specialisation holds no member variable and many functions are NO-OPs because if the alphabet has only
 * one valid value there is no state that can be changed.
 */
template <typename derived_type, typename char_t>
class alphabet_base<derived_type, 1ul, char_t>
{
protected:
    /*!\name Member types
     * \{
     */
    //!\copybrief seqan3::alphabet_base::char_type
    using char_type = std::conditional_t<std::same_as<char_t, void>, char, char_t>;
    //!\copybrief seqan3::alphabet_base::rank_type
    using rank_type = bool;
    //!\}

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr alphabet_base()                                   noexcept = default; //!< Defaulted.
    constexpr alphabet_base(alphabet_base const &)              noexcept = default; //!< Defaulted.
    constexpr alphabet_base(alphabet_base &&)                   noexcept = default; //!< Defaulted.
    constexpr alphabet_base & operator=(alphabet_base const &)  noexcept = default; //!< Defaulted.
    constexpr alphabet_base & operator=(alphabet_base &&)       noexcept = default; //!< Defaulted.
    ~alphabet_base()                                            noexcept = default; //!< Defaulted.
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\copybrief seqan3::alphabet_base::to_char
    constexpr char_type to_char() const noexcept
    //!\cond
        requires (!std::same_as<char_t, void>)
    //!\endcond
    {
        return derived_type::char_value;
    }

    //!\copybrief seqan3::alphabet_base::to_rank
    constexpr rank_type to_rank() const noexcept
    {
        return 0;
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\copybrief seqan3::alphabet_base::assign_char
    constexpr derived_type & assign_char(char_type const) noexcept
    //!\cond
        requires (!std::same_as<char_t, void>)
    //!\endcond
    {
        return static_cast<derived_type &>(*this);
    }

    //!\copybrief seqan3::alphabet_base::assign_rank
    constexpr derived_type & assign_rank(rank_type const) noexcept
    {
        return static_cast<derived_type &>(*this);
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr bool alphabet_size = 1;

    //!\name Comparison operators
    //!\{

    //!\brief Letters are always equal.
    friend constexpr bool operator==(derived_type const, derived_type const) noexcept
    {
        return true;
    }

    //!\brief Letters are never unequal.
    friend constexpr bool operator!=(derived_type const, derived_type const) noexcept
    {
        return false;
    }

    //!\brief One letter cannot be smaller than another.
    friend constexpr bool operator<(derived_type const,  derived_type const)  noexcept
    {
        return false;
    }

    //!\brief One letter cannot be bigger than another.
    friend constexpr bool operator>(derived_type const,  derived_type const)  noexcept
    {
        return false;
    }

    //!\brief Letters are always equal.
    friend constexpr bool operator<=(derived_type const, derived_type const) noexcept
    {
        return true;
    }

    //!\brief Letters are always equal.
    friend constexpr bool operator>=(derived_type const, derived_type const) noexcept
    {
        return true;
    }
    //!\}

private:
#if SEQAN3_WORKAROUND_GCC_87113
    bool _bug_workaround{};
#endif
};

} // namespace seqan3
