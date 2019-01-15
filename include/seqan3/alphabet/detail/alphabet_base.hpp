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
#include <seqan3/core/detail/int_types.hpp>
#include <seqan3/std/concepts>

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
    //TODO: this should be > 1 and there should be a different specialisation for size 1
    static_assert(size > 0, "It does not make sense to use the base class for alphabets of size < 1.");

    /*!\name Member types
     * \{
     */
    //!\brief The type of the alphabet when converted to char (e.g. via to_char()).
    using char_type = char_t;
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
    /*!\brief Assign from a character.
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
    constexpr derived_type & assign_char(std::conditional_t<std::Same<char_type, void>, char, char_type> const c) noexcept
    //!\cond
        requires !std::Same<char_type, void>
    //!\endcond
    {
        using index_t = std::make_unsigned_t<char_type>;
        rank = derived_type::char_to_rank[static_cast<index_t>(c)];
        return static_cast<derived_type &>(*this);
    }

    /*!\brief Assign from a numeric value.
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

private:
    //!\brief The value of the alphabet letter is stored as the rank.
    rank_type rank;

    //!\brief Expose the derived_type back to the derived_type as a dependent type to allow some obscure template
    //!       witchery.
    using derived_t = derived_type;

    //!\brief Enable usage of the previous typedef.
    friend derived_type;
};

} // namespace seqan3
