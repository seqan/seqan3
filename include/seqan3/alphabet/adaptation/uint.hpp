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

/*!\file alphabet/adaptation/uint.hpp
 * \ingroup adaptation
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides alphabet adaptations for standard uint types.
 * \details
 * This file provides function and metafunction overloads so that the following types
 * fulfil the seqan3::alphabet_concept:
 *   * `uint`
 *   * `uint16_t`
 *   * `uint32_t`
 *
 * You will likely not use these interfaces directly, they are, however, very helpful
 * for conversions between other alphabets and between other alphabets and characters.
 *
 * \attention
 * Note that `uint64_t` is absent from the list, because of there is no corresponding
 * character type.
 *
 * \attention
 * Please be aware that if you also include alphabet/concept.hpp, you need to do so **after**
 * including this file, not before.
 */

#pragma once

#include <limits>

#include <meta/meta.hpp>

#include <seqan3/alphabet/concept_fwd.hpp>

//!\cond
// NOTE(h-2): 'tis some nice code witchery that gives an error if the include order is wrong.
namespace
{

template <std::nullptr_t t>
concept bool alphabet_concept = true;

}

namespace seqan3
{

static_assert(alphabet_concept<nullptr>, "You must include alphabet/concept.hpp AFTER including alphabet/adaptation/uint.hpp, not before!");

}
//!\endcond

namespace seqan3::detail
{
//!\addtogroup adaptation
//!\{

//!\brief The list of types that are defined as uint adaptations.
using uint_adaptations            = meta::list<uint8_t,
                                               uint16_t,
                                               uint32_t>;
//!\brief The corresponding list of rank types
using uint_adaptations_char_types = meta::list<char,
                                               char16_t,
                                               char32_t>;

//!\brief Metafunction that indicates whether a type is a uint alphabet adaptation.
template <typename type>
struct is_uint_adaptation :
    public std::conditional_t<meta::in<uint_adaptations, type>::value, std::true_type, std::false_type>
{};

//!\brief Shortcut for seqan3::detail::is_uint_adaptation.
template <typename type>
constexpr bool is_uint_adaptation_v = is_uint_adaptation<type>::value;

//!\}
} // namespace seqan3::detail

namespace seqan3
{

//!\addtogroup adaptation
//!\{

// ------------------------------------------------------------------
// type metafunctions
// ------------------------------------------------------------------

//!\brief Type metafunction that shall return the type of an alphabet in uint representation. [uint specialization]
//!\tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
//!\sa seqan3::underlying_char
template <typename uint_type>
    requires detail::is_uint_adaptation_v<uint_type>
struct underlying_char<uint_type>
{
    //!\brief The character type of the same size as `uint_type`.
    using type = meta::at<detail::uint_adaptations_char_types,
                          meta::find_index<detail::uint_adaptations, uint_type>>;
};

//!\brief Type metafunction that shall return the type of an alphabet in rank representation. [uint specialization]
//!\tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
//!\sa seqan3::underlying_rank
template <typename uint_type>
    requires detail::is_uint_adaptation_v<uint_type>
struct underlying_rank<uint_type>
{
    //!\brief The same as `uint_type`.
    using type = uint_type;
};

// ------------------------------------------------------------------
// value metafunctions
// ------------------------------------------------------------------

//!\brief Value metafunction that shall return the size of an alphabet. [uint specialization]
//!\tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
//!\sa seqan3::alphabet_size
template <typename uint_type>
    requires detail::is_uint_adaptation_v<uint_type>
struct alphabet_size<uint_type>
{
    //!\brief Smallest unsigned integral type that can hold value;
    using type = detail::min_viable_uint_t<(std::size_t)std::numeric_limits<uint_type>::max() + 1 -
                                            std::numeric_limits<uint_type>::lowest()>;
    //!\brief The alphabet's size.
    static constexpr type value =
        (type)std::numeric_limits<uint_type>::max() + 1 - std::numeric_limits<uint_type>::lowest();
};

// ------------------------------------------------------------------
// free functions
// ------------------------------------------------------------------

/*!\name Free function wrappers for the uint alphabet adaptation
 * \brief For `uint8_t`, `uint16_t` and `uint32_t` do conversion to/from char types.
 * \ingroup adaptation
 * \{
 */
/*!\brief Converting uint to char casts to an character type of same size.
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param c The alphabet letter that you wish to convert to char.
 * \returns The letter's value in the alphabet's rank type (usually `uint`).
 */
template <typename uint_type>
constexpr underlying_char_t<uint_type> to_char(uint_type const c)
    requires detail::is_uint_adaptation_v<uint_type>
{
    return c;
}

/*!\brief Converting uint to rank is a no-op (it will just return the value you pass in).
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param c The alphabet letter that you wish to convert to rank.
 * \returns The letter's value in the alphabet's rank type (usually a `uint*_t`).
 */
template <typename uint_type>
constexpr underlying_rank_t<uint_type> to_rank(uint_type const c)
    requires detail::is_uint_adaptation_v<uint_type>
{
    return c;
}

/*!\brief Assign from a character type vie implicit or explicit cast.
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param c The alphabet letter that you wish to assign to.
 * \param in The `uint` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename uint_type>
constexpr uint_type & assign_char(uint_type & c, underlying_char_t<uint_type> const in)
    requires detail::is_uint_adaptation_v<uint_type>
{
    return c = in;
}

/*!\brief Assign from a character type vie implicit or explicit cast.
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param c An alphabet letter temporary.
 * \param in The `uint` value you wish to assign.
 * \returns The assignment result as a temporary.
 * \details
 * Use this e.g. to newly create alphabet letters from uint:
 * ```cpp
 * auto l = assign_char(dna5{}, 'G');  // l is of type dna5
 * ```
 */
template <typename uint_type>
constexpr uint_type && assign_char(uint_type && c, underlying_char_t<uint_type> const in)
    requires detail::is_uint_adaptation_v<std::remove_reference_t<uint_type>>
{
    return std::move(c = in);
}

/*!\brief The same as calling `=` on the uint.
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param c The alphabet letter that you wish to assign to.
 * \param in The `rank` value you wish to assign.
 * \returns A reference to the alphabet letter you passed in.
 */
template <typename uint_type>
constexpr uint_type & assign_rank(uint_type & c, underlying_rank_t<uint_type> const in)
    requires detail::is_uint_adaptation_v<uint_type>
{
    return c = in;
}

/*!\brief The same as calling `=` on the uint.
 * \tparam uint_type One of `uint8_t`, `uint16_t` or `uint32_t`.
 * \param c An alphabet letter temporary.
 * \param in The `rank` value you wish to assign.
 * \returns The assignment result as a temporary.
 * \details
 * Use this e.g. to newly create alphabet letters from rank:
 * ```cpp
 * auto l = assign_rank(dna5{}, 1);  // l is of type dna5 and == dna5::C
 * ```
 */
template <typename uint_type>
constexpr uint_type && assign_rank(uint_type && c, underlying_rank_t<uint_type> const in)
    requires detail::is_uint_adaptation_v<std::remove_reference_t<uint_type>>
{
    return std::move(c = in);
}
//!\}
//!\}

} // namespace seqan3

// concept tests in alphabet/adaptation.hpp because of include order constraints
