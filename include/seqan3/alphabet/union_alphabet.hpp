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

#include "nucleotide/dna5.hpp"
#include "alphabet.hpp"

/*! The union alphabet
 * \ingroup alphabet
 */

namespace seqan3
{

namespace detail
{
template <typename tuple_t, size_t ...I>
constexpr auto alphabet_objects(std::index_sequence<I...>)
{
    using variant_t = std::variant<std::tuple_element_t<I, tuple_t>...>;
    using list_t = std::initializer_list<variant_t>;
    return list_t
    {
        variant_t
        {
            std::in_place_index_t<I>{},
            std::tuple_element_t<I, tuple_t>{}
        }...
    };
}

template <typename tuple_t, size_t ...I>
constexpr auto alphabet_sizes(std::index_sequence<I...>)
{
    using integral_t = uint8_t;
    using list_t = std::initializer_list<integral_t>;
    return list_t{std::tuple_element_t<I, tuple_t>::value_size...};
}

template <size_t size>
constexpr auto at_least_type() {
    if constexpr(size <= 1)
    {
        return bool{};
    } else if constexpr(size <= 255u)
    {
        return uint8_t{};
    } else if constexpr(size <= 65535u)
    {
        return uint16_t{};
    } else if constexpr(size <= 4294967295u)
    {
        return uint32_t{};
    } else
    {
        return uint64_t{};
    }
}

} // namespace seqan::detail

/*! The union alphabet
 *
 * The union alphabet represents the union of two or more alphabets (e.g. the four letter DNA alphabet +
 * the gap alphabet). The alphabet may be brace initialized from the static letter members (see above).
 * Note that you cannot assign regular characters, but additional functions for this are available.
 *
 *     union_alphabet<dna4, gap> my_letter{};
 *     // doesn't work:
 *     // union_alphabet<dna4, gap> my_letter{'A'};
 *
 *     union_alphabet<dna4, gap>.from_char('C'); // <- this does!
 *
 *     union_alphabet<dna4, gap>.from_char('-'); // gap character
 *     if (my_letter.to_char() == 'A')
 *        std::cout << "yeah\n"; // "yeah";
 */

template <typename first_alphabet_type, typename ...alphabet_types>
    requires alphabet_concept<first_alphabet_type> && (alphabet_concept<alphabet_types> && ...)
struct union_alphabet
{
    /* types */
public:
    //! The size of the alphabet, i.e. the number of different values it can take.
    static constexpr size_t value_size = (alphabet_types::value_size + ... + first_alphabet_type::value_size);
    //! the type of the alphabet when converted to char (e.g. via @link to_char @endlink)
    using char_type = typename first_alphabet_type::char_type;
    //! the type of the alphabet when represented as a number (e.g. via @link to_integral @endlink)
    using integral_type = decltype(detail::at_least_type<value_size>());

protected:
    using tuple_t = std::tuple<first_alphabet_type, alphabet_types...>;
    using index_sequence_t = decltype(std::make_index_sequence<std::tuple_size_v<tuple_t>>{});
    static constexpr auto sizes = detail::alphabet_sizes<tuple_t>(index_sequence_t{});
    static constexpr auto alphabets = detail::alphabet_objects<tuple_t>(index_sequence_t{});

public:
    //! internal value
    integral_type value;

    //! ability to cast to @link char_type @endlink **explicitly**.
    explicit constexpr operator char_type() const
    {
        return to_char();
    }

    //! return the letter as a character of @link char_type @endlink.
    constexpr char_type to_char() const
    {
        auto sizes_it = sizes.begin();
        auto sizes_it_end = sizes.end();
        auto alphabets_it = alphabets.begin();
        auto size_start = 0;
        auto size_end = 0;

        std::optional<char_type> result;

        while(sizes_it != sizes_it_end && !result.has_value()) {
            size_end += *sizes_it;
            std::visit([&](auto alphabet){
                if(value < size_end) {
                    result = alphabet.from_integral(value - size_start).to_char();
                }
            }, *alphabets_it);

            size_start = size_end;
            sizes_it++;
            alphabets_it++;
        }

        return result.value_or(static_cast<char_type>(0));
    }

    //! return the letter's numeric value or rank in the alphabet
    constexpr integral_type to_integral() const
    {
        return value;
    }

    //! assign from a character
    constexpr union_alphabet & from_char(char_type const c)
    {
        auto sizes_it = sizes.begin();
        auto sizes_it_end = sizes.end();
        auto alphabets_it = alphabets.begin();
        auto size_start = 0;

        bool found = false;

        while(sizes_it != sizes_it_end && !found) {
            std::visit([&](auto alphabet){
                alphabet.from_char(c);
                if(alphabet.to_char() == c) {
                    value = size_start + alphabet.to_integral();
                    found = true;
                }
            }, *alphabets_it);

            size_start += *sizes_it;
            sizes_it++;
            alphabets_it++;
        }

        if(!found) {
            value = 0;
        }

        return *this;
    }

    //! assign from a numeric value
    constexpr union_alphabet & from_integral(integral_type const i)
    {
        assert(i < value_size);
        value = i;
        return *this;
    }

    //! @name comparison operators
    //!@{
    constexpr bool operator==(union_alphabet const & rhs) const
    {
        return value == rhs.value;
    }

    constexpr bool operator!=(union_alphabet const & rhs) const
    {
        return value != rhs.value;
    }

    constexpr bool operator<(union_alphabet const & rhs) const
    {
        return value < rhs.value;
    }

    constexpr bool operator>(union_alphabet const & rhs) const
    {
        return value > rhs.value;
    }

    constexpr bool operator<=(union_alphabet const & rhs) const
    {
        return value <= rhs.value;
    }

    constexpr bool operator>=(union_alphabet const & rhs) const
    {
        return value >= rhs.value;
    }
    //!@}
};

static_assert(detail::internal_alphabet_concept<union_alphabet<dna5, dna5>>);
static_assert(alphabet_concept<union_alphabet<dna5, dna5>>);

} // namespace seqan3
