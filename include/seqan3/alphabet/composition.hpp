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
#include <seqan3/core/detail/int_types.hpp>

/*!\file alphabet/composition.hpp
 * \ingroup alphabet
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>

 * \brief Contains alphabet_composition.
 */

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

    //!\brief An index sequence up to the number of contained letters.
    static constexpr auto _positions = std::make_index_sequence<sizeof...(alphabet_types)>{};
    //!\endcond


    /*!\name Read functions
     * \{
     */
    /*!\brief Return the letter combinations numeric value or rank in the alphabet composition.
     * \par Complexity
     * Linear in the number of alpahabets.
     */
    constexpr integral_type to_integral() const
    {
        return to_integral_impl(_positions);
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    /*!\brief Assign from a numeric value.
     * \par Complexity
     * Linear in the number of alpahabets.
     * \par Exceptions
     * Asserts that the parameter is smaller than value_size [only in debug mode].
     */
    constexpr derived_type & from_integral(integral_type const i)
    {
        assert(i < value_size);
        from_integral_impl<0>(i);
        return static_cast<derived_type &>(*this);
    }
    //!\}

    //\brief Assign to the first letter that matches in the composition.
    //TODO this is not passed on to descendents, why?
//     template <typename type>
//     constexpr derived_type & operator=(type && val)
//         requires meta::in<meta::list<first_alphabet_type, alphabet_types...>, type>::value
//     {
//         std::get<type>(*this) = std::forward<type>(val);
//         return static_cast<derived_type &>(*this);
//     }


private:
    //!\brief declared private to prevent direct use of the CRTP base
    alphabet_composition() = default;
    //!\brief declared private to prevent direct use of the CRTP base
    constexpr alphabet_composition(alphabet_composition const &) = default;
    //!\brief declared private to prevent direct use of the CRTP base
    constexpr alphabet_composition(alphabet_composition &&) = default;
    //!\brief declared private to prevent direct use of the CRTP base
    constexpr alphabet_composition & operator =(alphabet_composition const &) = default;
    //!\brief declared private to prevent direct use of the CRTP base
    constexpr alphabet_composition & operator =(alphabet_composition &&) = default;
    //!\brief declared private to prevent direct use of the CRTP base
    ~alphabet_composition() = default;

    //!\brief befriend the derived type so that it can instantiate
    //!\sa https://isocpp.org/blog/2017/04/quick-q-prevent-user-from-derive-from-incorrect-crtp-base
    friend derived_type;

    //!\brief Implementation of to_integral().
    template <std::size_t ...idx>
    constexpr integral_type to_integral_impl(std::index_sequence<idx...> const &) const
    {
        if constexpr (sizeof...(idx) > 0)
        {
            return seqan3::to_integral(std::get<0>(*this)) +
                   ((seqan3::to_integral(std::get<idx + 1>(*this)) * _cummulative_alph_sizes[idx]) + ...);
        }
        else
        {
            return seqan3::to_integral(std::get<0>(*this));
        }
    }

    //!\brief Implementation of from_integral().
    template <std::size_t j>
    constexpr void
    from_integral_impl(integral_type const i)
    {
        if constexpr (j == 0)
        {
            seqan3::from_integral(std::get<j>(*this),
                                i % alphabet_size_v<meta::at_c<meta::list<first_alphabet_type, alphabet_types...>, j>>);
        } else
        {
            seqan3::from_integral(std::get<j>(*this),
                          (i / _cummulative_alph_sizes[j-1])
                          % alphabet_size_v<meta::at_c<meta::list<first_alphabet_type, alphabet_types...>, j>>);
        }

        if constexpr (j < sizeof...(alphabet_types))
            from_integral_impl<j+1>(i);
    }
};

} // namespace seqan3
