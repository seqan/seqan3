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

#include <iostream>
#include <string>
#include <utility>

#include <seqan3/alphabet/alphabet.hpp>
#include <seqan3/alphabet/detail/pod_tuple.hpp>

/*!\file alphabet/alphabet_composition.hpp
 * \ingroup alphabet
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>

 * \brief Contains alphabet_composition.
 */

namespace seqan3
{

template <typename first_alphabet_type, typename ...alphabet_types>
      requires alphabet_concept<first_alphabet_type> && (alphabet_concept<alphabet_types> && ...)
struct alphabet_composition;

} // namespace seqan3

namespace seqan3::detail
{

template <uint64_t value>
using min_viable_uint_t = std::conditional_t<value < 255ull,        uint8_t,
                          std::conditional_t<value < 65535ull,      uint16_t,
                          std::conditional_t<value < 4294967295ull, uint32_t, uint64_t>>>;

} // namespace seqan3::detail

namespace seqan3::detail::alphabet_composition
{

template <std::size_t i, typename alphabet_one, typename ...alphabet_types>
constexpr typename alphabet_composition<alphabet_one, alphabet_types...>::integral_type
cummulative_alph_size()
{
    static_assert(sizeof...(alphabet_types) >= i);

    if constexpr (i == 0)
       return alphabet_size_v<alphabet_one>;
    else
       return alphabet_size_v<alphabet_one> * cummulative_alph_size<i-1, alphabet_types...>();

}

// try to figure this one out :D
template <std::size_t ...idx,
          typename ...alphabet_types>
constexpr typename alphabet_composition<alphabet_types...>::integral_type
to_integral_impl(std::index_sequence<idx...> const &,
                 seqan3::alphabet_composition<alphabet_types...> const & comp)
{
    if constexpr (sizeof...(idx) > 0)
    {
        return to_integral(std::get<0>(comp)) +
               ((to_integral(std::get<idx + 1>(comp)) * cummulative_alph_size<idx, alphabet_types...>()) +  ...);
    }
    else
    {
        return to_integral(std::get<0>(comp));
    }
}

template <std::size_t j, typename ...alphabet_types>
constexpr void
from_integral_impl(seqan3::alphabet_composition<alphabet_types...> & comp,
                   typename seqan3::alphabet_composition<alphabet_types...>::integral_type const i)
{
    if constexpr (j == 0)
    {
        from_integral(std::get<j>(comp),
                      i % alphabet_size_v<get_ith_type_t<j, alphabet_types...>>);
    } else
    {
        from_integral(std::get<j>(comp),
                      (i / cummulative_alph_size<j-1, alphabet_types...>())
                       % alphabet_size_v<get_ith_type_t<j, alphabet_types...>>);
    }

    if constexpr (j + 1 < sizeof...(alphabet_types))
        from_integral_impl<j+1>(comp, i);
}

} // namespace seqan3::detail::alphabet_composition


namespace seqan3
{

/*!\brief The basis of alphabets that contain multiple (different) letters at one position.
 * \ingroup alphabet
 * \tparam first_alphabet_type Type of the first letter; must satisfy alphabet_concept.
 * \tparam alphabet_types Types of further letters (up to 4); must satisfy alphabet_concept.
 *
 * This data structure provides the basis of a combined alphabet, where the different
 * alphabet letters exist independently, similar to a tuple. In fact this class
 * provides a tuple-like interface with `std::get<0>(t)` and objects can be brace-initialized
 * with the individual members.
 *
 * TODO example
 *
 *
 * An alphabet_composition itself does not satisfy the alphabet_concept. If you desire this
 * behaviour (usually you do), inherit from it and add `from_char` and `to_char` member functions.
 * See quality_composition or mask_composition for more details.
 */

template <typename first_alphabet_type, typename ...alphabet_types>
      requires alphabet_concept<first_alphabet_type> && (alphabet_concept<alphabet_types> && ...)
struct alphabet_composition : public detail::pod_tuple<first_alphabet_type, alphabet_types...>
{
    //!\brief The type of value_size and `alphabet_size_v<alphabet_composition<...>>`
    using integral_type = detail::min_viable_uint_t<alphabet_size_v<first_alphabet_type> *
                                                    (alphabet_size_v<alphabet_types> * ...)>;

    //!\brief The product of the sizes of the individual alphabets.
    static constexpr integral_type value_size{alphabet_size_v<first_alphabet_type> *
                                              (alphabet_size_v<alphabet_types> * ...)};

    //!\privatesection
    static constexpr auto _positions = std::make_index_sequence<sizeof...(alphabet_types)>{};
    //!\publicsection

    //! TODO
    constexpr integral_type to_integral() const
    {

        return detail::alphabet_composition::to_integral_impl(_positions, *this);
    }

    //! TODO
    constexpr alphabet_composition & from_integral(integral_type const i)
    {
        detail::alphabet_composition::from_integral_impl<0>(*this, i);
        return *this;
    }
};

}

// } // namespace seqan3
