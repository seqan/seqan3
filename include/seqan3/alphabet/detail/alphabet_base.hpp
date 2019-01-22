// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::alphabet_base.
 */

#pragma once

#include <seqan3/alphabet/detail/member_exposure.hpp>
#include <seqan3/alphabet/exception.hpp>
#include <seqan3/core/detail/int_types.hpp>
#include <seqan3/core/detail/reflection.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/type_traits>

namespace seqan3
{

/*!\brief A CRTP-base that makes defining a custom alphabet easier.
 * \ingroup alphabet
 * \tparam derived_type The CRTP parameter type.
 * \tparam size         The size of the alphabet.
 * \tparam char_t       The character type of the alphabet (set this to `void` when defining just a
 *                      seqan3::semi_alphabet_concept).
 *
 * \details
 *
 * You can use this class to define your own alphabet, but types are not required to be based on it to model
 * seqan3::alphabet_concept, it is purely a way to avoid code duplication.
 *
 * The base class represents the alphabet value as the rank and
 * automatically deduces the rank type from the size, it further defines all required member functions and types; the
 * derived type needs to define only the following two tables as static member variables:
 *
 *   * `static std::array<char_type, value_size> constexpr rank_to_char` that defines for every possible rank value
 *     the corresponding char value.
 *   * `static std::array<rank_type, 256> constexpr char_to_rank` that defines for every possible character value the
 *     corresponding rank value.
 *
 * ### Example
 *
 * This creates an alphabet called `ab` which has size two and the two letters 'A' and 'B':
 * \snippet test/snippet/alphabet/detail/alphabet_base.cpp example
 *
 */
template <typename derived_type, size_t size, typename char_t = char>
class alphabet_base
{
public:
    static_assert(size > 1, "It does not make sense to use the base class for alphabets of size <= 1.");

    /*!\name Member types
     * \{
     */
    //!\brief The type of the alphabet when converted to char (e.g. via to_char()).
    using char_type = char_t;
    //!\brief Equal to char_type in all relevant situations (needed only to make semi alphabet definitions legal).
    using char_type_ = std::conditional_t<std::Same<char_type, void>, char, char_type>;
    //!\brief The type of the alphabet when represented as a number (e.g. via to_rank()).
    using rank_type = detail::min_viable_uint_t<size - 1>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr alphabet_base() noexcept : rank{} {}
    constexpr alphabet_base(alphabet_base const &) = default;
    constexpr alphabet_base(alphabet_base &&) = default;
    constexpr alphabet_base & operator=(alphabet_base const &) = default;
    constexpr alphabet_base & operator=(alphabet_base &&) = default;
    ~alphabet_base() = default;
    //!\}

    /*!\name Read functions
     * \{
     */
    /*!\brief Return the letter as a character of char_type.
     *
     * \details
     *
     * Satisfies the seqan3::alphabet_concept::to_char() requirement via the seqan3::to_char() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr char_type to_char() const noexcept
    //!\cond
        requires !std::Same<char_type, void>
    //!\endcond
    {
        return derived_type::rank_to_char[rank];
    }

    /*!\brief Return the letter's numeric value (rank in the alphabet).
     *
     * \details
     *
     * Satisfies the seqan3::semi_alphabet_concept::to_rank() requirement via the to_rank() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
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
     * Satisfies the seqan3::alphabet_concept::assign_char() requirement via the seqan3::assign_char() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr derived_type & assign_char(char_type_ const c) noexcept
    //!\cond
        requires !std::Same<char_type, void>
    //!\endcond
    {
        using index_t = std::make_unsigned_t<char_type>;
        rank = derived_type::char_to_rank[static_cast<index_t>(c)];
        return static_cast<derived_type &>(*this);
    }

    /*!\brief Assign from a character, throw on invalid characters.
     * \param c The character to be assigned.
     * \throws invalid_char_assignment If seqan3::char_is_valid
     *
     * \details
     *
     * Satisfies the seqan3::alphabet_concept::assign_char_strict() requirement via the seqan3::assign_char_strict()
     * wrapper.
     *
     * ### Complexity
     *
     * Constant; but slightly slower than #assign_char, because it performs checks.
     */
    derived_type & assign_char_strict(char_type_ const c)
    //!\cond
        requires !std::Same<char_type, void>
    //!\endcond
    {
        if (!derived_type::char_is_valid(c))
            throw invalid_char_assignment{detail::get_display_name_v<derived_type>.string(), c};

        using seqan3::assign_char;
        return assign_char(static_cast<derived_type &>(*this), c);
    }

    /*!\brief Assign from a numeric value.
     * \param c The rank to be assigned.
     *
     * \details
     *
     * Satisfies the seqan3::semi_alphabet_concept::assign_rank() requirement via the seqan3::assign_rank() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    constexpr derived_type & assign_rank(rank_type const c) noexcept
    {
        assert(static_cast<size_t>(c) < static_cast<size_t>(value_size));
        rank = c;
        return static_cast<derived_type &>(*this);
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static detail::min_viable_uint_t<size> constexpr value_size = size;

    //!\name Comparison operators
    //!\{
    friend constexpr bool operator==(derived_type const & lhs, derived_type const & rhs) noexcept
    {
        using seqan3::to_rank;
        return to_rank(lhs) == to_rank(rhs);
    }

    friend constexpr bool operator!=(derived_type const & lhs, derived_type const & rhs) noexcept
    {
        using seqan3::to_rank;
        return to_rank(lhs) != to_rank(rhs);
    }

    friend constexpr bool operator<(derived_type const & lhs, derived_type const & rhs) noexcept
    {
        using seqan3::to_rank;
        return to_rank(lhs) < to_rank(rhs);
    }

    friend constexpr bool operator>(derived_type const & lhs, derived_type const & rhs) noexcept
    {
        using seqan3::to_rank;
        return to_rank(lhs) > to_rank(rhs);
    }

    friend constexpr bool operator<=(derived_type const & lhs, derived_type const & rhs) noexcept
    {
        using seqan3::to_rank;
        return to_rank(lhs) <= to_rank(rhs);
    }

    friend constexpr bool operator>=(derived_type const & lhs, derived_type const & rhs) noexcept
    {
        using seqan3::to_rank;
        return to_rank(lhs) >= to_rank(rhs);
    }
    //!\}

    /*!\brief Validate whether a character value has a one-to-one mapping to an alphabet value.
     *
     * \details
     *
     * Models the seqan3::semi_alphabet_concept::char_is_valid_for() requirement via the seqan3::char_is_valid_for()
     * wrapper.
     *
     * Default implementation: True for all character values that are reproduced by #to_char() after being assigned
     * to the alphabet.
     *
     * \par Complexity
     *
     * Constant.
     *
     * \par Exceptions
     *
     * Guaranteed not to throw.
     */
    static constexpr bool char_is_valid(char_type_ const c) noexcept
    //!\cond
        requires !std::Same<char_type, void>
    //!\endcond
    {
        using seqan3::to_char;
        using seqan3::assign_char;
        return to_char(assign_char(derived_type{}, c)) == c;
    }

private:
    //!\brief The value of the alphabet letter is stored as the rank.
    rank_type rank;
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
public:
    /*!\name Member types
     * \{
     */
    //!\copybrief seqan3::alphabet_base::char_type
    using char_type = char_t;
    //!\copybrief seqan3::alphabet_base::char_type_
    using char_type_ = std::conditional_t<std::Same<char_type, void>, char, char_type>;
    //!\copybrief seqan3::alphabet_base::rank_type
    using rank_type = bool;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr alphabet_base() noexcept = default;
    constexpr alphabet_base(alphabet_base const &) = default;
    constexpr alphabet_base(alphabet_base &&) = default;
    constexpr alphabet_base & operator=(alphabet_base const &) = default;
    constexpr alphabet_base & operator=(alphabet_base &&) = default;
    ~alphabet_base() = default;
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\copybrief seqan3::alphabet_base::to_char
    constexpr char_type to_char() const noexcept
    //!\cond
        requires !std::Same<char_type, void>
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
    constexpr derived_type & assign_char(char_type_ const) noexcept
    //!\cond
        requires !std::Same<char_type, void>
    //!\endcond
    {
        return static_cast<derived_type &>(*this);
    }

    //!\copybrief seqan3::alphabet_base::assign_char_strict
    derived_type & assign_char_strict(char_type_ const c)
    //!\cond
        requires !std::Same<char_type, void>
    //!\endcond
    {
        if (!derived_type::char_is_valid(c))
            throw invalid_char_assignment{detail::get_display_name_v<derived_type>.string(), c};

        using seqan3::assign_char;
        return assign_char(static_cast<derived_type &>(*this), c);
    }

    //!\copybrief seqan3::alphabet_base::assign_rank
    constexpr derived_type & assign_rank(rank_type const) noexcept
    {
        return static_cast<derived_type &>(*this);
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr bool value_size = 1;

    //!\name Comparison operators
    //!\{
    friend constexpr bool operator==(derived_type const &, derived_type const &) noexcept
    {
        return true;
    }

    friend constexpr bool operator!=(derived_type const &, derived_type const &) noexcept
    {
        return false;
    }

    friend constexpr bool operator<(derived_type const &,  derived_type const &)  noexcept
    {
        return false;
    }

    friend constexpr bool operator>(derived_type const &,  derived_type const &)  noexcept
    {
        return false;
    }

    friend constexpr bool operator<=(derived_type const &, derived_type const &) noexcept
    {
        return true;
    }

    friend constexpr bool operator>=(derived_type const &, derived_type const &) noexcept
    {
        return true;
    }
    //!\}

    //!\copybrief seqan3::alphabet_base::char_is_valid
    static constexpr bool char_is_valid(char_type_ const c) noexcept
    //!\cond
        requires !std::Same<char_type, void>
    //!\endcond
    {
        return derived_type::char_value == c;
    }

private:
    //!\cond
    bool _bug_workaround; // See GCC Bug-Report: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=87113
    //!\endcond
};

} // namespace seqan3
