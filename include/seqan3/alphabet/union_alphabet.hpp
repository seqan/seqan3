// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#pragma once

#include <tuple>
#include <variant>
#include <utility>
#include <optional>
#include <cassert>
#include <algorithm>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/int_types.hpp>

/*! The union alphabet
 * \ingroup alphabet
 */

namespace seqan3::detail
{

//!\privatesection

/*!\brief Returns an array which contains the prefix sum over all alphabet_types::value_size's.
 * \relates seqan3::union_alphabet
 *
 * ```cpp
 * constexpr auto prefix_sum = alphabet_prefix_sum_sizes<dna4, gap, dna5>();
 * assert(prefix_sum.size() == 4);
 * assert(prefix_sum[0] == 0);
 * assert(prefix_sum[1] == 4);
 * assert(prefix_sum[2] == 5);
 * assert(prefix_sum[3] == 10);
 * ```
 */
template <typename ...alphabet_types>
constexpr auto alphabet_prefix_sum_sizes()
{
    constexpr size_t value_size = ((size_t)alphabet_types::value_size + ... + 0);
    using rank_t = min_viable_uint_t<value_size>;

    constexpr size_t N = sizeof...(alphabet_types) + 1;
    using array_t = std::array<rank_t, N>;

    array_t prefix_sum{0, alphabet_types::value_size...};
    for (auto i = 1u; i < N; ++i)
        prefix_sum[i] += prefix_sum[i-1];

    return prefix_sum;
}
//!\publicsection

} // namespace seqan3::detail

namespace seqan3::detail::union_alphabet
{

//!\privatesection

/*!\brief Returns an fixed sized map at compile time where the key is the rank
 * of alphabet_t and the value is the corresponding char of that rank.
 * \relates seqan3::detail::union_alphabet::value_to_char_table
 *
 * ```cpp
 * constexpr auto table1 = value_to_char_table_I<5, char>(dna4{});
 * assert(table1.size() == 5);
 * assert(table1[0] == 'A');
 * assert(table1[1] == 'C');
 * assert(table1[2] == 'G');
 * assert(table1[3] == 'T');
 * assert(table1[4] == '\0');
 * ```
 */
template <size_t max_value_size, typename char_t, typename alphabet_t>
constexpr auto value_to_char_table_I(alphabet_t alphabet)
{
    using array_t = std::array<char_t, max_value_size>;
    array_t value_to_char_{};

    for (auto i = 0; i < alphabet_t::value_size; ++i)
        value_to_char_[i] = alphabet.assign_rank(i).to_char();

    return value_to_char_;
}

/*!\brief Returns an map at compile time where the key is the rank of the union
 * of all alphabets and the value is the corresponding char of that rank and
 * alphabet.
 * \relates seqan3::union_alphabet
 *
 * ```cpp
 * constexpr auto value_to_char = value_to_char_table<char, dna4, gap, dna5>();
 * assert(value_to_char.size() == 10);
 * assert(value_to_char[0] == 'A');
 * assert(value_to_char[1] == 'C');
 * assert(value_to_char[2] == 'G');
 * assert(value_to_char[3] == 'T');
 * assert(value_to_char[4] == '-');
 * assert(value_to_char[5] == 'A');
 * assert(value_to_char[6] == 'C');
 * // and so on
 * ```
 */
template <typename char_t, typename ...alphabet_types>
constexpr auto value_to_char_table()
{
    constexpr auto table_size = (alphabet_types::value_size + ... + 0);
    constexpr auto value_sizes = std::array<size_t, table_size>{alphabet_types::value_size...};
    constexpr auto max_value_size
        = std::max({static_cast<size_t>(0), static_cast<size_t>(alphabet_types::value_size)...});

    using array_t = std::array<char_t, table_size>;
    using array_inner_t = std::array<char_t, max_value_size>;
    using array_array_t = std::array<array_inner_t, table_size>;

    constexpr auto array_array = array_array_t
    {
        value_to_char_table_I<max_value_size, char_t>(alphabet_types{})...
    };

    array_t value_to_char{};
    for (auto i = 0u, value = 0u; i < table_size; ++i)
        for (auto k = 0u; k < value_sizes[i]; ++k, ++value)
            value_to_char[value] = array_array[i][k];

    return value_to_char;
}

/*!\brief Returns an map at compile time where the key is the char of one of the
 * alphabets and the value is the corresponding rank over all alphabets (by
 * conflict will default to the first).
 * \relates seqan3::union_alphabet
 *
 * ```cpp
 * constexpr auto char_to_value = char_to_value_table<char, dna4, gap, dna5>();
 * assert(char_to_value.size() == 256);
 * assert(char_to_value['A'] == 0);
 * assert(char_to_value['C'] == 1);
 * assert(char_to_value['G'] == 2);
 * assert(char_to_value['T'] == 3);
 * assert(char_to_value['-'] == 4);
 * assert(char_to_value['A'] == 0);
 * assert(char_to_value['C'] == 1);
 * assert(char_to_value['G'] == 2);
 * assert(char_to_value['T'] == 3);
 * assert(char_to_value['N'] == 9);
 * assert(char_to_value['*'] == 0); // every other character defaults to 0
 * ```
 */
template <typename char_t, typename ...alphabet_types>
constexpr auto char_to_value_table()
{
    constexpr size_t value_size = ((size_t)alphabet_types::value_size + ... + 0);
    using rank_t = min_viable_uint_t<value_size>;

    constexpr auto table_size = 1 << (sizeof(char_t) * 8);
    constexpr auto value_to_char = value_to_char_table<char_t, alphabet_types...>();

    using array_t = std::array<rank_t, table_size>;

    array_t char_to_value{};
    for (auto i = 0u; i < value_to_char.size(); ++i)
    {
        auto & old_entry = char_to_value[value_to_char[i]];
        auto is_new_entry = value_to_char[0] != value_to_char[i] && old_entry == 0;
        if (is_new_entry)
            old_entry = static_cast<rank_t>(i);
    }
    return char_to_value;
}
//!\publicsection

} // namespace seqan3::detail::union_alphabet

namespace seqan3
{

/*!\brief An union_alphabet that merges different regular alphabets as a single alphabet.
 * \ingroup alphabet
 * \tparam first_alphabet_type Type of the first letter, e.g. dna4; must satisfy seqan3::alphabet_concept.
 * \tparam alphabet_types Types of further letters; must satisfy seqan3::alphabet_concept.
 *
 * The union alphabet represents the union of two or more alphabets (e.g. the
 * four letter DNA alphabet + the gap alphabet). Note that you cannot assign
 * regular characters, but additional functions for this are available.
 *
 * ```cpp
 *     union_alphabet<dna4, gap> my_letter{};
 *     union_alphabet<dna4, gap> converted_letter{dna4::C};
 *     // doesn't work:
 *     // union_alphabet<dna4, gap> my_letter{'A'};
 *
 *     union_alphabet<dna4, gap>.assign_char('C'); // <- this does!
 *     union_alphabet<dna4, gap>.assign_char('-'); // gap character
 *     union_alphabet<dna4, gap>.assign_char('K'); // unknown characters map to the default/unknown
 *                                               // character of the first alphabet type (i.e. A of dna4)
 *     if (my_letter.to_char() == 'A')
 *        std::cout << "yeah\n"; // "yeah";
 * ```
 *
 * The union alphabet can also be assigned or constructed directly from one of
 * the base alphabets.
 *
 * ```cpp
 * using alphabet_t = union_alphabet<dna4, dna5, gap>;
 * using variant_t = alphabet_t::variant_type;
 *
 * constexpr alphabet_t letter0{dna4::A};
 * constexpr alphabet_t letter1 = {dna4::C};
 * constexpr alphabet_t letter2 = static_cast<variant_t>(dna4::G);
 * constexpr alphabet_t letter3{dna4::T}; // letter3 = dna4::T; does not work
 *
 * assert(letter0.to_rank() == 0);
 * assert(letter1.to_rank() == 1);
 * assert(letter2.to_rank() == 2);
 * assert(letter3.to_rank() == 3);
 *
 * alphabet_t letter4{dna5::A}; // letter4 = dna5::A; does not work
 * alphabet_t letter5 = {dna5::C};
 * alphabet_t letter6 = static_cast<variant_t>(dna5::G);
 *
 * assert(letter4.to_rank() == 4);
 * assert(letter5.to_rank() == 5);
 * assert(letter6.to_rank() == 6);
 * ```
 */
template <typename first_alphabet_type, typename ...alphabet_types>
    requires alphabet_concept<first_alphabet_type> && (alphabet_concept<alphabet_types> && ...)
class union_alphabet
{
public:
    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr size_t value_size = (alphabet_types::value_size + ... + first_alphabet_type::value_size);

    //!\brief The type of the alphabet when converted to char (e.g. via \link to_char \endlink)
    using char_type = typename first_alphabet_type::char_type;

    //!\brief The type of the alphabet when represented as a number (e.g. via \link to_rank \endlink)
    using rank_type = detail::min_viable_uint_t<value_size>;

    /*!\brief The type used to assign a value from one of the base alphabets
     * during copy construction or copy assignment.
     *
     * For example used by
     * seqan3::union_alphabet::union_alphabet(const variant_type & alphabet) or
     * seqan3::union_alphabet::operator=(const variant_type & alphabet)
     */
    using variant_type = std::variant<first_alphabet_type, alphabet_types...>;

    //!\name Default constructors
    //!\{
    constexpr union_alphabet() = default;
    constexpr union_alphabet(union_alphabet const &) = default;
    constexpr union_alphabet(union_alphabet &&) = default;
    //!\}

    //!\name Default assignment operators
    //!\{
    constexpr union_alphabet & operator= (union_alphabet const &) = default;
    constexpr union_alphabet & operator= (union_alphabet &&) = default;
    //!\}

    //!\name Conversion constructors
    //\{
    explicit constexpr union_alphabet(rank_type && value) :
        _value{value}
    {}
    explicit constexpr union_alphabet(rank_type const & value) :
        _value{value}
    {}

    /*!\brief Construction via a value of the base alphabets
     *
     * ```cpp
     *     union_alphabet<dna4, gap> letter1{dna4::C}; // or
     *     union_alphabet<dna4, gap> letter2 = gap::GAP;
     * ```
     */
    constexpr union_alphabet(variant_type const & alphabet) :
        _value{from_base_(alphabet)}
    {}
    //!\}

    //!\name Conversion assignment operators
    //!\{
    /*!\brief Assignment via a value of the base alphabets
     *
     * ```cpp
     *     union_alphabet<dna4, gap> letter1{};
     *     letter1 = gap::GAP;
     * ```
     */
    constexpr union_alphabet & operator= (variant_type const & alphabet)
    {
        _value = from_base_(alphabet);
        return *this;
    }
    //!\}

    //!\brief Ability to cast to \link char_type \endlink **explicitly**.
    explicit constexpr operator char_type() const
    {
        return to_char();
    }

    /*!\name Read functions
     * \{
     */
    //!\brief Return the letter as a character of char_type.
    constexpr char_type to_char() const
    {
        return value_to_char[_value];
    }

    //!\brief Return the letter's numeric value or rank in the alphabet.
    constexpr rank_type to_rank() const
    {
        return _value;
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\brief Assign from a character.
    constexpr union_alphabet & assign_char(char_type const c)
    {
        _value = char_to_value[c];
        return *this;
    }

    //!\brief Assign from a numeric value.
    constexpr union_alphabet & assign_rank(rank_type const i)
    {
        assert(i < value_size);
        _value = i;
        return *this;
    }
    //!\}

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(union_alphabet const & rhs) const
    {
        return _value == rhs._value;
    }

    constexpr bool operator!=(union_alphabet const & rhs) const
    {
        return _value != rhs._value;
    }

    constexpr bool operator<(union_alphabet const & rhs) const
    {
        return _value < rhs._value;
    }

    constexpr bool operator>(union_alphabet const & rhs) const
    {
        return _value > rhs._value;
    }

    constexpr bool operator<=(union_alphabet const & rhs) const
    {
        return _value <= rhs._value;
    }

    constexpr bool operator>=(union_alphabet const & rhs) const
    {
        return _value >= rhs._value;
    }
    //!\}

    //!\privatesection
    //!\brief The data member.
    rank_type _value;
    //!\publicsection

protected:
    //!\privatesection

    //!\brief Converts an object of one of the base alphabets into the internal representation
    static constexpr rank_type from_base_(variant_type const & alphabet_v)
    {
        return std::visit([&](auto && alphabet) -> rank_type
        {
            return prefix_sum_sizes[alphabet_v.index()] + static_cast<rank_type>(alphabet.to_rank());
        }, alphabet_v);
    }

    //!\brief Compile-time generated lookup table which contains the prefix sum up to the position of each alphabet
    //!\hideinitializer
    //!\sa alphabet_prefix_sum_sizes
    static constexpr auto prefix_sum_sizes =
        detail::alphabet_prefix_sum_sizes<first_alphabet_type, alphabet_types...>();

    //!\brief Compile-time generated lookup table which maps the rank to char
    //!\hideinitializer
    //!\sa value_to_char_table
    static constexpr auto value_to_char =
        detail::union_alphabet::value_to_char_table<char_type, first_alphabet_type, alphabet_types...>();

    //!\brief Compile-time generated lookup table which maps the char to rank
    //!\hideinitializer
    //!\sa char_to_value_table
    static constexpr auto char_to_value =
        detail::union_alphabet::char_to_value_table<char_type, first_alphabet_type, alphabet_types...>();
};

} // namespace seqan3

#ifndef NDEBUG
#include "nucleotide/dna5.hpp"
static_assert(seqan3::alphabet_concept<seqan3::union_alphabet<seqan3::dna5, seqan3::dna5>>);
#endif
