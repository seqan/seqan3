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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// Author: Marcel Ehrhardt <marcel.ehrhardt@fu-berlin.de>
// Author: David Heller <david.heller@fu-berlin.de>
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

} // namespace seqan3::detail

namespace seqan3::detail::union_alphabet
{

template <size_t max_value_size, typename char_t, typename alphabet_t>
constexpr auto value_to_char_table_I(alphabet_t alphabet)
{
    using array_t = std::array<char_t, max_value_size>;
    array_t value_to_char_{};

    for (auto i = 0; i < alphabet_t::value_size; ++i)
        value_to_char_[i] = alphabet.assign_rank(i).to_char();

    return value_to_char_;
}

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

} // namespace seqan3::detail::union_alphabet

namespace seqan3
{

/*!\brief An union_alphabet that merges different regular alphabets as a single alphabet.
 * \ingroup alphabet
 * \tparam first_alphabet_type Type of the first letter, e.g. dna4; must satisfy seqan3::alphabet_concept.
 * \tparam alphabet_types Types of further letters; must satisfy seqan3::alphabet_concept.
 *
 * The union alphabet represents the union of two or more alphabets (e.g. the four letter DNA alphabet +
 * the gap alphabet). The alphabet may be brace initialized from the static letter members (see above).
 * Note that you cannot assign regular characters, but additional functions for this are available.
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
 */
template <typename first_alphabet_type, typename ...alphabet_types>
    requires alphabet_concept<first_alphabet_type> && (alphabet_concept<alphabet_types> && ...)
class union_alphabet
{
public:
    /* types */
    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr size_t value_size = (alphabet_types::value_size + ... + first_alphabet_type::value_size);
    //!\brief the type of the alphabet when converted to char (e.g. via \link to_char \endlink)
    using char_type = typename first_alphabet_type::char_type;
    //!\brief the type of the alphabet when represented as a number (e.g. via \link to_rank \endlink)
    using rank_type = detail::min_viable_uint_t<value_size>;
    /*!\brief the type used to assign a value from one of the base alphabets
     * during copy construction or copy assignment.
     *
     * (i.e. via
     * \link union_alphabet(const variant_type & alphabet) \endlink or
     * \link operator=(const variant_type & alphabet) \endlink )
     */
    using variant_type = std::variant<first_alphabet_type, alphabet_types...>;

    //!\name default constructors
    //!\{
    constexpr union_alphabet() = default;
    constexpr union_alphabet(union_alphabet const &) = default;
    constexpr union_alphabet(union_alphabet &&) = default;
    //!\}

    //!\name default assignment operators
    //!\{
    constexpr union_alphabet & operator= (union_alphabet const &) = default;
    constexpr union_alphabet & operator= (union_alphabet &&) = default;
    //!\}

    //!\name default assignment operators
    //!\{
    /*explicit*/ constexpr union_alphabet(rank_type && value)
        : _value{value}
    {}
    /*explicit*/ constexpr union_alphabet(rank_type const & value)
        : _value{value}
    {}
    //!\}

    /*!\brief allow construction via a value of the base alphabets
     *
     * ```cpp
     *     union_alphabet<dna4, gap> letter1{dna4::C}; // or
     *     union_alphabet<dna4, gap> letter2 = gap::GAP;
     * ```
     */
    constexpr union_alphabet(variant_type const & alphabet)
        : _value{from_base_(alphabet)}
    {}

    /*!\brief allow assignment via a value of the base alphabets
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

    //! internal value
    rank_type _value;

    //! ability to cast to \link char_type \endlink **explicitly**.
    explicit constexpr operator char_type() const
    {
        return to_char();
    }

    //! return the letter as a character of \link char_type \endlink.
    constexpr char_type to_char() const
    {
        return value_to_char[_value];
    }

    //! return the letter's numeric value or rank in the alphabet
    constexpr rank_type to_rank() const
    {
        return _value;
    }

    //! assign from a character
    constexpr union_alphabet & assign_char(char_type const c)
    {
        _value = char_to_value[c];
        return *this;
    }

    //! assign from a numeric value
    constexpr union_alphabet & assign_rank(rank_type const i)
    {
        assert(i < value_size);
        _value = i;
        return *this;
    }

    //! @name comparison operators
    //!@{
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
    //!@}

protected:
    //! \privatesection
    // conversion tables

    static constexpr rank_type from_base_(variant_type const & alphabet_v)
    {
        return std::visit([&](auto && alphabet) -> rank_type
        {
            return prefix_sum_sizes[alphabet_v.index()] + static_cast<rank_type>(alphabet.to_rank());
        }, alphabet_v);
    }

    static constexpr auto prefix_sum_sizes
        = detail::alphabet_prefix_sum_sizes<first_alphabet_type, alphabet_types...>();

    static constexpr auto value_to_char
        = detail::union_alphabet::value_to_char_table<char_type, first_alphabet_type, alphabet_types...>();

    static constexpr auto char_to_value
        = detail::union_alphabet::char_to_value_table<char_type, first_alphabet_type, alphabet_types...>();
};

} // namespace seqan3

#ifndef NDEBUG
#include "nucleotide/dna5.hpp"
static_assert(std::is_pod_v<seqan3::union_alphabet<seqan3::dna5, seqan3::dna5>>);
static_assert(std::is_trivial_v<seqan3::union_alphabet<seqan3::dna5, seqan3::dna5>>);
static_assert(std::is_trivially_copyable_v<seqan3::union_alphabet<seqan3::dna5, seqan3::dna5>>);
static_assert(std::is_standard_layout_v<seqan3::union_alphabet<seqan3::dna5, seqan3::dna5>>);
static_assert(seqan3::alphabet_concept<seqan3::union_alphabet<seqan3::dna5, seqan3::dna5>>);
#endif
