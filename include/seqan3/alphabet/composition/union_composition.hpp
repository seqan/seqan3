// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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

/*!\brief A composition that merges different regular alphabets into a single alphabet.
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
public:
    /*!\brief Returns true if alphabet_t is one of the given alphabet types.
     * \tparam alphabet_t The type to check.
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
    static constexpr auto value_size =
        detail::min_viable_uint_v<
            (static_cast<size_t>(alphabet_size_v<alphabet_types>) + ... + static_cast<size_t>(0))
        >;

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
     * \tparam alphabet_t One of the base alphabet types.
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
     * \tparam I The index of the i-th base alphabet.
     * \tparam alphabet_t The i-th given base alphabet type.
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
     * \tparam alphabet_t One of the base alphabet types.
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
        using index_t = std::make_unsigned_t<char_type>;
        _value = char_to_value[static_cast<index_t>(c)];
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

    /*!\brief Compile-time generated lookup table which contains the partial
     * sum up to the position of each alphabet.
     *
     * An array which contains the prefix sum over all
     * alphabet_types::value_size's.
     *
     * ```cpp
     * constexpr std::array partial_sum = union_composition<dna4, gap, dna5>::partial_sum_sizes; // not working; is protected
     * assert(partial_sum.size() == 4);
     * assert(partial_sum[0] == 0);
     * assert(partial_sum[1] == 4);
     * assert(partial_sum[2] == 5);
     * assert(partial_sum[3] == 10);
     * ```
     */
    static constexpr std::array partial_sum_sizes = []() constexpr
    {
        constexpr size_t N = sizeof...(alphabet_types) + 1;

        std::array<rank_type, N> partial_sum{0, alphabet_size_v<alphabet_types>...};
        for (size_t i = 1u; i < N; ++i)
            partial_sum[i] += partial_sum[i-1];

        return partial_sum;
    }();

    /*!\brief Compile-time generated lookup table which maps the rank to char.
     *
     * A map generated at compile time where the key is the rank of the union
     * of all alphabets and the value is the corresponding char of that rank
     * and alphabet.
     *
     * ```cpp
     * constexpr std::array value_to_char = union_composition<char, dna4, gap, dna5>::value_to_char; // not working; is protected
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
    static constexpr std::array<char_type, value_size> value_to_char = []() constexpr
    {
        // Explicitly writing assign_rank_to_char within assign_value_to_char
        // causes this bug (g++-7 and g++-8):
        // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=84684
        auto assign_rank_to_char = [](auto alphabet, size_t rank) constexpr
        {
            return seqan3::to_char(seqan3::assign_rank(alphabet, rank));
        };

        auto assign_value_to_char = [assign_rank_to_char] (auto alphabet, auto & value_to_char, auto & value) constexpr
        {
            using alphabet_t = std::decay_t<decltype(alphabet)>;
            for (size_t i = 0u; i < alphabet_size_v<alphabet_t>; ++i, ++value)
                value_to_char[value] = assign_rank_to_char(alphabet, i);
        };

        unsigned value = 0u;
        std::array<char_type, value_size> value_to_char{};

        // initializer lists guarantee sequencing;
        // the following expression behaves as:
        // for(auto alphabet: alphabet_types)
        //    assign_value_to_char(alphabet, value_to_char, value);
        ((assign_value_to_char(alphabet_types{}, value_to_char, value)),...);

        return value_to_char;
    }();

    /*!\brief Compile-time generated lookup table which maps the char to rank.
     *
     * An map generated at compile time where the key is the char of one of the
     * alphabets and the value is the corresponding rank over all alphabets (by
     * conflict will default to the first).
     *
     * ```cpp
     * constexpr std::array char_to_value = char_to_value_table<char, dna4, gap, dna5>();
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
    static constexpr std::array char_to_value = []() constexpr
    {
        constexpr size_t table_size = 1 << (sizeof(char_type) * 8);

        std::array<rank_type, table_size> char_to_value{};
        for (size_t i = 0u; i < value_to_char.size(); ++i)
        {
            using index_t = std::make_unsigned_t<char_type>;
            rank_type & old_entry = char_to_value[static_cast<index_t>(value_to_char[i])];
            bool is_new_entry = value_to_char[0] != value_to_char[i] && old_entry == 0;
            if (is_new_entry)
                old_entry = static_cast<rank_type>(i);
        }
        return char_to_value;
    }();

    //!\brief Converts an object of one of the given alphabets into the internal representation.
    //!\tparam index The position of `alphabet_t` in the template pack `alphabet_types`.
    //!\tparam alphabet_t One of the base alphabet types.
    //!\param alphabet The value of a base alphabet that should be assigned.
    template <size_t index, typename alphabet_t>
    //!\cond
        requires has_alternative<alphabet_t>()
    //!\endcond
    static constexpr rank_type rank_by_index_(alphabet_t const & alphabet) noexcept
    {
        return partial_sum_sizes[index] + static_cast<rank_type>(seqan3::to_rank(alphabet));
    }

    //!\brief Converts an object of one of the given alphabets into the internal representation.
    //!\details Finds the index of alphabet_t in the given types.
    //!\tparam alphabet_t One of the base alphabet types.
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
