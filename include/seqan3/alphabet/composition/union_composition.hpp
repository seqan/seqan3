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

/*!\file
 * \ingroup composition
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \author David Heller <david.heller AT fu-berlin.de>
 * \brief Contains seqan3::union_composition.
 */

#pragma once

#include <array>
#include <utility>
#include <cassert>
#include <algorithm>

#include <meta/meta.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/int_types.hpp>

namespace seqan3
{

/*!\brief A composition that merges different regular alphabets as a single alphabet.
 * \ingroup composition
 * \tparam ...alphabet_types Types of further letters; must satisfy seqan3::alphabet_concept, e.g. dna4.
 * \implements seqan3::alphabet_concept
 *
 * The union alphabet represents the union of two or more alphabets (e.g. the
 * four letter DNA alphabet + the gap alphabet). Note that you cannot assign
 * regular characters, but additional functions for this are available.
 *
 * This class has a similar behavior as std::variant.
 *
 * ```cpp
 *     union_composition<dna4, gap> my_letter{};
 *     union_composition<dna4, gap> converted_letter{dna4::C};
 *     // doesn't work:
 *     // union_composition<dna4, gap> my_letter{'A'};
 *
 *     union_composition<dna4, gap>{}.assign_char('C'); // <- this does!
 *     union_composition<dna4, gap>{}.assign_char('-'); // gap character
 *     union_composition<dna4, gap>{}.assign_char('K'); // unknown characters map to the default/unknown
 *                                                   // character of the first alphabet type (i.e. A of dna4)
 *     if (my_letter.to_char() == 'A')
 *        std::cout << "yeah\n"; // "yeah";
 * ```
 *
 * The union alphabet can also be constructed directly from one of the base
 * alphabets.
 *
 * ```cpp
 * using alphabet_t = union_composition<dna4, dna5, gap>;
 *
 * constexpr alphabet_t letter0{dna4::A};
 * constexpr alphabet_t letter1 = dna4::C;
 * constexpr alphabet_t letter2 = {dna4::G};
 * constexpr alphabet_t letter3 = static_cast<alphabet_t>(dna4::T);
 *
 * assert(letter0.to_rank() == 0);
 * assert(letter1.to_rank() == 1);
 * assert(letter2.to_rank() == 2);
 * assert(letter3.to_rank() == 3);
 * ```
 *
 * Or can be assigned by one of the base alphabets.
 *
 * ```cpp
 * using alphabet_t = union_composition<dna4, dna5, gap>;
 *
 * alphabet_t letter;
 *
 * letter = dna5::A;
 * assert(letter.to_rank() == 4);
 *
 * letter = {dna5::C};
 * assert(letter.to_rank() == 5);
 *
 * letter = static_cast<alphabet_t>(dna5::G);
 * assert(letter.to_rank() == 6);
 * ```
 */
template <typename ...alphabet_types>
//!\cond
    requires (alphabet_concept<alphabet_types> && ... && (sizeof...(alphabet_types) >= 1))
//!\endcond
class union_composition
{
protected:
    //!\privatesection

    /*!\brief Returns the sum over all alphabet_size_v<alphabets_t>'s.
     * \tparam ...alphabets_t The types must satisfy seqan3::alphabet_concept.
     * \hideinitializer
     *
     * ```cpp
     * constexpr auto sum = sum_of_alphabet_sizes_v<dna4, gap, dna5>;
     * assert(sum == 10);
     * ```
     */
    template <typename ...alphabets_t>
    //!\cond
        requires (alphabet_concept<alphabets_t> && ...)
    //!\endcond
    static constexpr auto sum_of_alphabet_sizes_v =
        (static_cast<size_t>(alphabet_size_v<alphabets_t>) + ... + static_cast<size_t>(0));

    //!\publicsection
public:
    /*!\brief Returns true if alphabet_t is one of the given alphabet types.
     * \tparam alphabet_t The type to check
     *
     * ```cpp
     * using union_t = union_composition<dna4, gap>;
     *
     * static_assert(union_t::has_alternative<dna4>(), "should be true");
     * static_assert(union_t::has_alternative<gap>(), "should be true");
     * static_assert(!union_t::has_alternative<dna5>(), "should be false");
     * ```
     */
    template <typename alphabet_t>
    static constexpr bool has_alternative() noexcept
    {
        return meta::in<meta::list<alphabet_types...>, alphabet_t>::value;
    }

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr size_t value_size = sum_of_alphabet_sizes_v<alphabet_types...>;

    //!\brief The type of the alphabet when converted to char (e.g. via \link to_char \endlink).
    using char_type = underlying_char_t<meta::front<meta::list<alphabet_types...>>>;

    //!\brief The type of the alphabet when represented as a number (e.g. via \link to_rank \endlink).
    using rank_type = detail::min_viable_uint_t<value_size>;

    /*!\name Default constructors
     * \{
     */
    constexpr union_composition() = default;
    constexpr union_composition(union_composition const &) = default;
    constexpr union_composition(union_composition &&) = default;
    //!\}

    /*!\name Default assignment operators
     * \{
     */
    constexpr union_composition & operator= (union_composition const &) = default;
    constexpr union_composition & operator= (union_composition &&) = default;
    //!\}

    /*!\name Conversion constructors
     * \{
     */

    /*!\brief Construction via a value of the base alphabets.
     * \tparam alphabet_t One of the base alphabet types
     * \param alphabet The value of a base alphabet that should be assigned.
     *
     * ```cpp
     *     union_composition<dna4, gap> letter1{dna4::C}; // or
     *     union_composition<dna4, gap> letter2 = gap::GAP;
     * ```
     */
    template <typename alphabet_t>
    //!\cond
        requires has_alternative<alphabet_t>()
    //!\endcond
    constexpr union_composition(alphabet_t const & alphabet) noexcept :
        _value{rank_by_type_(alphabet)}
    {}

    /*!\brief Construction via a value of reoccurring alphabets.
     * \tparam I The index of the i-th base alphabet
     * \tparam alphabet_t The i-th given base alphabet type
     * \param alphabet The value of a base alphabet that should be assigned.
     *
     * ```cpp
     * using alphabet_t = union_composition<dna4, dna4>;
     *
     * constexpr alphabet_t letter0{std::in_place_index_t<0>{}, dna4::A};
     * constexpr alphabet_t letter4{std::in_place_index_t<1>{}, dna4::A};
     *
     * EXPECT_EQ(letter0.to_rank(), 0);
     * EXPECT_EQ(letter4.to_rank(), 4);
     * ```
     */
    template <size_t I, typename alphabet_t>
    //!\cond
        requires has_alternative<alphabet_t>()
    //!\endcond
    constexpr union_composition(std::in_place_index_t<I>, alphabet_t const & alphabet) noexcept :
        _value{rank_by_index_<I>(alphabet)}
    {}
    //!\}

    /*!\name Conversion assignment operators
     * \{
     */

    /*!\brief Assignment via a value of the base alphabets.
     * \tparam alphabet_t One of the base alphabet types
     * \param alphabet The value of a base alphabet that should be assigned.
     *
     * ```cpp
     *     union_composition<dna4, gap> letter1{};
     *     letter1 = gap::GAP;
     * ```
     */
    template <typename alphabet_t>
    //!\cond
        requires has_alternative<alphabet_t>()
    //!\endcond
    constexpr union_composition & operator= (alphabet_t const & alphabet) noexcept
    {
        _value = rank_by_type_(alphabet);
        return *this;
    }
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\brief Return the letter as a character of char_type.
    constexpr char_type to_char() const noexcept
    {
        return value_to_char[_value];
    }

    //!\brief Return the letter's numeric value or rank in the alphabet.
    constexpr rank_type to_rank() const noexcept
    {
        return _value;
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\brief Assign from a character.
    //!\param c The character to assign to.
    constexpr union_composition & assign_char(char_type const c) noexcept
    {
        _value = char_to_value[c];
        return *this;
    }

    //!\brief Assign from a numeric value.
    //!\param i The rank of a character to assign to.
    constexpr union_composition & assign_rank(rank_type const i) /*noexcept*/
    {
        // TODO(marehr): mark function noexcept if assert is replaced
        // https://github.com/seqan/seqan3/issues/85
        assert(i < value_size);
        _value = i;
        return *this;
    }
    //!\}

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(union_composition const & rhs) const noexcept
    {
        return _value == rhs._value;
    }

    constexpr bool operator!=(union_composition const & rhs) const noexcept
    {
        return _value != rhs._value;
    }

    constexpr bool operator<(union_composition const & rhs) const noexcept
    {
        return _value < rhs._value;
    }

    constexpr bool operator>(union_composition const & rhs) const noexcept
    {
        return _value > rhs._value;
    }

    constexpr bool operator<=(union_composition const & rhs) const noexcept
    {
        return _value <= rhs._value;
    }

    constexpr bool operator>=(union_composition const & rhs) const noexcept
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

    /*!\brief Returns the max over all alphabet_size_v<alphabets_t>'s.
     * \tparam ...alphabets_t The types must satisfy seqan3::alphabet_concept.
     *
     * ```cpp
     * constexpr auto max = max_of_alphabet_sizes_v<dna4, gap, dna5>;
     * assert(max == 5);
     * ```
     */
    template <typename ...alphabets_t>
    //!\cond
        requires (alphabet_concept<alphabets_t> && ...)
    //!\endcond
    static constexpr auto max_of_alphabet_sizes_v =
        std::max({static_cast<size_t>(0), static_cast<size_t>(alphabet_size_v<alphabets_t>)...});

    /*!\brief Returns an array which contains the prefix sum over all alphabet_types::value_size's.
     * \tparam ...alphabet_types The types must satisfy seqan3::alphabet_concept.
     *
     * ```cpp
     * constexpr auto partial_sum = partial_sum_of_alphabet_sizes<dna4, gap, dna5>();
     * assert(partial_sum.size() == 4);
     * assert(partial_sum[0] == 0);
     * assert(partial_sum[1] == 4);
     * assert(partial_sum[2] == 5);
     * assert(partial_sum[3] == 10);
     * ```
     */
    template <typename ...alphabets_t>
    //!\cond
        requires (alphabet_concept<alphabets_t> && ...)
    //!\endcond
    static constexpr auto partial_sum_of_alphabet_sizes() noexcept
    {
        constexpr auto value_size = sum_of_alphabet_sizes_v<alphabets_t...>;
        using rank_t = detail::min_viable_uint_t<value_size>;

        constexpr size_t N = sizeof...(alphabets_t) + 1;
        using array_t = std::array<rank_t, N>;

        array_t partial_sum{0, alphabet_size_v<alphabets_t>...};
        for (auto i = 1u; i < N; ++i)
            partial_sum[i] = static_cast<rank_t>(partial_sum[i] + partial_sum[i-1]);

        return partial_sum;
    }

    /*!\brief Returns an fixed sized map at compile time where the key is the rank
     * of alphabet_t and the value is the corresponding char of that rank.
     * \sa value_to_char_table
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
    template <size_t max_value_size, typename alphabet_t>
    static constexpr auto value_to_char_table_I(alphabet_t alphabet) noexcept
    {
        using array_t = std::array<char_type, max_value_size>;
        array_t value_to_char_{};

        for (char_type i = 0; i < alphabet_size_v<alphabet_t>; ++i)
            value_to_char_[i] = alphabet.assign_rank(i).to_char();

        return value_to_char_;
    }

    /*!\brief Returns an map at compile time where the key is the rank of the union
     * of all alphabets and the value is the corresponding char of that rank and
     * alphabet.
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
    static constexpr auto value_to_char_table() noexcept
    {
        constexpr auto table_size = sum_of_alphabet_sizes_v<alphabet_types...>;
        constexpr auto value_sizes = std::array<size_t, table_size>
        {
            alphabet_size_v<alphabet_types>...
        };
        constexpr auto max_value_size = max_of_alphabet_sizes_v<alphabet_types...>;

        using array_t = std::array<char_type, table_size>;
        using array_inner_t = std::array<char_type, max_value_size>;
        using array_array_t = std::array<array_inner_t, table_size>;

        constexpr auto array_array = array_array_t
        {
            value_to_char_table_I<max_value_size>(alphabet_types{})...
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
    static constexpr auto char_to_value_table() noexcept
    {
        constexpr auto value_size = sum_of_alphabet_sizes_v<alphabet_types...>;
        using rank_t = detail::min_viable_uint_t<value_size>;

        constexpr auto table_size = 1 << (sizeof(char_type) * 8);
        constexpr auto value_to_char = value_to_char_table();

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

    //!\brief Compile-time generated lookup table which contains the prefix sum up to the position of each alphabet.
    //!\sa partial_sum_of_alphabet_sizes
    static constexpr auto partial_sum_sizes = partial_sum_of_alphabet_sizes<alphabet_types...>();

    //!\brief Compile-time generated lookup table which maps the rank to char.
    //!\sa value_to_char_table
    static constexpr auto value_to_char = value_to_char_table();

    //!\brief Compile-time generated lookup table which maps the char to rank.
    //!\sa char_to_value_table
    static constexpr auto char_to_value = char_to_value_table();

    //!\brief Converts an object of one of the given alphabets into the internal representation.
    //!\tparam index The position of `alphabet_t` in the template pack `alphabet_types`
    //!\tparam alphabet_t One of the base alphabet types
    //!\param alphabet The value of a base alphabet that should be assigned.
    template <size_t index, typename alphabet_t>
    //!\cond
        requires has_alternative<alphabet_t>()
    //!\endcond
    static constexpr rank_type rank_by_index_(alphabet_t const & alphabet) noexcept
    {
        return static_cast<rank_type>(partial_sum_sizes[index] + alphabet.to_rank());
    }

    //!\brief Converts an object of one of the given alphabets into the internal representation.
    //!\details Finds the index of alphabet_t in the given types.
    //!\tparam alphabet_t One of the base alphabet types
    //!\param alphabet The value of a base alphabet that should be assigned.
    template <typename alphabet_t>
    //!\cond
        requires has_alternative<alphabet_t>()
    //!\endcond
    static constexpr rank_type rank_by_type_(alphabet_t const & alphabet) noexcept
    {
        constexpr size_t index = meta::find_index<meta::list<alphabet_types...>, alphabet_t>::value;
        return rank_by_index_<index>(alphabet);
    }
};

} // namespace seqan3
