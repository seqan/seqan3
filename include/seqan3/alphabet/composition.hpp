// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================
// Author: Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
// ============================================================================

#pragma once

#include <cassert>
#include <utility>

#include <meta/meta.hpp>

#include <seqan3/alphabet/alphabet.hpp>
#include <seqan3/alphabet/quality/concept.hpp>
#include <seqan3/core/pod_tuple.hpp>

/*!\file alphabet/composition.hpp
 * \ingroup alphabet
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>

 * \brief Contains alphabet_composition.
 */

namespace seqan3
{

template <typename derived_type,
          typename first_alphabet_type,
          typename ...alphabet_types>
      requires alphabet_concept<first_alphabet_type> && (alphabet_concept<alphabet_types> && ...)
struct alphabet_composition;

} // namespace seqan3

namespace seqan3::detail
{

// TODO this needs to go to core somewhere
template <uint64_t value>
using min_viable_uint_t = std::conditional_t<value < 255ull,        uint8_t,
                          std::conditional_t<value < 65535ull,      uint16_t,
                          std::conditional_t<value < 4294967295ull, uint32_t, uint64_t>>>;

} // namespace seqan3::detail

namespace seqan3::detail::alphabet_composition
{

template <std::size_t ...idx,
          typename derived_type,
          typename ...alphabet_types>
constexpr typename alphabet_composition<derived_type, alphabet_types...>::integral_type
to_integral_impl(std::index_sequence<idx...> const &,
                 seqan3::alphabet_composition<derived_type, alphabet_types...> const & comp)
{
    if constexpr (sizeof...(idx) > 0)
    {
        return to_integral(std::get<0>(comp)) +
               ((to_integral(std::get<idx + 1>(comp)) * seqan3::alphabet_composition<derived_type, alphabet_types...>::_cummulative_alph_sizes[idx]) + ...);
    }
    else
    {
        return to_integral(std::get<0>(comp));
    }
}

template <std::size_t j,
          typename derived_type,
          typename ...alphabet_types>
constexpr void
from_integral_impl(seqan3::alphabet_composition<derived_type, alphabet_types...> & comp,
                   typename seqan3::alphabet_composition<derived_type, alphabet_types...>::integral_type const i)
{
    if constexpr (j == 0)
    {
        from_integral(std::get<j>(comp),
                      i % alphabet_size_v<meta::at_c<meta::list<alphabet_types...>, j>>);
    } else
    {
        from_integral(std::get<j>(comp),
                      (i / seqan3::alphabet_composition<derived_type, alphabet_types...>::_cummulative_alph_sizes[j-1])
                       % alphabet_size_v<meta::at_c<meta::list<alphabet_types...>, j>>);
    }

    if constexpr (j + 1 < sizeof...(alphabet_types))
        from_integral_impl<j+1>(comp, i);
}

} // namespace seqan3::detail::alphabet_composition

namespace seqan3
{

/*!\brief The CRTP base of alphabets that contain multiple (different) letters at one position.
 * \ingroup alphabet
 * \tparam first_alphabet_type Type of the first letter; must satisfy alphabet_concept.
 * \tparam alphabet_types Types of further letters (up to 4); must satisfy alphabet_concept.
 *
 * This data structure is CRTP base clasee for combined alphabets, where the different
 * alphabet letters exist independently, similar to a tuple. In fact this class
 * provides a tuple-like interface with `std::get<0>(t)` and objects can be brace-initialized
 * with the individual members.
 *
 * \attention
 * This a "pure base class", you cannot instantiate it, you can only inherit from it.
 * Most likely you are interested in using one of it's descendents like quality_composition.
 *
 * \sa quality_composition
 * \sa mask_composition
 */

template <typename derived_type,
          typename first_alphabet_type,
          typename ...alphabet_types>
      requires alphabet_concept<first_alphabet_type> && (alphabet_concept<alphabet_types> && ...)
struct alphabet_composition :
    public pod_tuple<first_alphabet_type, alphabet_types...>
{
private:
    //!\brief declare private to prevent direct use of the CRTP base
    alphabet_composition() = default;
    //!\brief declare private to prevent direct use of the CRTP base
    constexpr alphabet_composition(alphabet_composition const &) = default;
    //!\brief declare private to prevent direct use of the CRTP base
    constexpr alphabet_composition(alphabet_composition &&) = default;
    //!\brief declare private to prevent direct use of the CRTP base
    constexpr alphabet_composition & operator =(alphabet_composition const &) = default;
    //!\brief declare private to prevent direct use of the CRTP base
    constexpr alphabet_composition & operator =(alphabet_composition &&) = default;
    //!\brief declare private to prevent direct use of the CRTP base
    ~alphabet_composition() = default;

    //!\brief befriend the derived type so that it can instantiate
    //!\sa https://isocpp.org/blog/2017/04/quick-q-prevent-user-from-derive-from-incorrect-crtp-base
    friend derived_type;

public:
    //!\brief The type of value_size and `alphabet_size_v<alphabet_composition<...>>`
    using integral_type = detail::min_viable_uint_t<alphabet_size_v<first_alphabet_type> *
                                                    (alphabet_size_v<alphabet_types> * ...)>;

    //!\brief The product of the sizes of the individual alphabets.
    static constexpr integral_type value_size{alphabet_size_v<first_alphabet_type> *
                                              (alphabet_size_v<alphabet_types> * ...)};

    //!\cond DEV
    //!\brief the cummulative alphabet size products (first, first*second, first*second*third...) are cached
    static constexpr std::array<integral_type, sizeof...(alphabet_types)+1> _cummulative_alph_sizes
    {
        [] () constexpr
        {
            std::array<integral_type, sizeof...(alphabet_types)+1> ret{};
            size_t count = 0;
            meta::for_each(meta::list<first_alphabet_type, alphabet_types...>{}, [&] (auto && alph) constexpr
            {
                ret[count] = alphabet_size_v<std::decay_t<decltype(alph)>> * (count > 0 ? ret[count - 1] : 1);
                ++count;
            });

            return std::move(ret);
        }()
    };
// this is more elegant, but harder to explain: 8-)
//     static constexpr std::array<integral_type, sizeof...(alphabet_types)+1> _cummulative_alph_sizes = meta::for_each(
//         meta::list<first_alphabet_type, alphabet_types...>{},
//         [ count = 0, ret = std::array<integral_type, sizeof...(alphabet_types)+1>{} ] (auto && alph) mutable constexpr
//     {
//         if constexpr(std::is_same_v<std::decay_t<decltype(alph)>, std::false_type>)
//         {
//             return std::move(ret);
//         }
//         else
//         {
//             ret[count] = alphabet_size_v<std::decay_t<decltype(alph)>> * (count > 0 ? ret[count - 1] : 1);
//             ++count;
//             return;
//         }
//     })(std::false_type{});

    //!\brief An index sequence up to the number of contained letters.
    static constexpr auto _positions = std::make_index_sequence<sizeof...(alphabet_types)>{};
    //!\endcond

    //! \brief Return the letter combinations numeric value or rank in the alphabet composition.
    constexpr integral_type to_integral() const
    {
        return detail::alphabet_composition::to_integral_impl(_positions, *this);
    }

    //! \brief Assign from a numeric value.
    constexpr derived_type & from_integral(integral_type const i)
    {
        assert(i < value_size);
        detail::alphabet_composition::from_integral_impl<0>(*this, i);
        return static_cast<derived_type &>(*this);
    }
};

} // namespace seqan3
