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

/*!\file
 * \brief Provides various metafunctions used by the range module.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <range/v3/utility/iterator_traits.hpp>
#include <range/v3/range_traits.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\addtogroup range
 * \{
 */

// ----------------------------------------------------------------------------

/*!\brief Whether a type has a `value_type` member and/or `ranges::value_type_t<t>` resolves [Value metafunction].
 * \tparam t The type to recurse on; must have `ranges::value_type_t<t>`
 * \hideinitializer
 */
template <typename t>
constexpr bool has_value_type_v = false;

//!\cond
template <typename t>
    requires requires (t) { typename ranges::value_type<std::decay_t<t>>::type; }
constexpr bool has_value_type_v<t> = true;
//!\endcond

// ----------------------------------------------------------------------------

/*!\brief Recursively determines the `value_type` on containers and/or iterators [Type metafunction].
 * \tparam t The type to recurse on; must have `ranges::value_type_t<t>`
 * \hideinitializer
 */
template <typename t>
//!\cond
    requires has_value_type_v<t>
//!\endcond
struct innermost_value_type
{
    //!\brief The forwarded type.
    using type = ranges::value_type_t<std::decay_t<t>>;
};

//!\cond
template <typename t>
    requires has_value_type_v<t> &&
             has_value_type_v<ranges::value_type_t<std::decay_t<t>>>
struct innermost_value_type<t>
{
    using type = typename innermost_value_type<ranges::value_type_t<std::decay_t<t>>>::type;
};
//!\endcond

// ----------------------------------------------------------------------------

//!\brief Shortcut for seqan3::innermost_value_type.
template <typename t>
//!\cond
    requires has_value_type_v<t>
//!\endcond
using innermost_value_type_t = typename innermost_value_type<t>::type;

// ----------------------------------------------------------------------------

/*!\brief Returns the number of times you can call `ranges::value_type_t` recursively on t [Value metafunction].
 * \tparam t The type to be queried; must resolve `ranges::value_type_t` at least once.
 * \hideinitializer
 */
template <typename t>
//!\cond
    requires has_value_type_v<t>
//!\endcond
constexpr size_t dimension_v = 1;

//!\cond
template <typename t>
    requires has_value_type_v<t> &&
             has_value_type_v<ranges::value_type_t<std::decay_t<t>>>
constexpr size_t dimension_v<t> = dimension_v<ranges::value_type_t<std::decay_t<t>>> + 1;
//!\endcond

// ----------------------------------------------------------------------------

/*!\interface seqan3::compatible_concept <>
 * \brief Two types are alike if their seqan3::dimension_v and their decayed seqan3::innermost_value_type_t are
 * the same.
 *
 * \details
 *
 * ```cpp
 * // these evaluate to true:
 * static_assert(seqan3::compatible_concept<std::string,              std::vector<char>>);
 * static_assert(seqan3::compatible_concept<std::vector<std::string>, std::vector<std::vector<char>>>);
 * ```
 */
//!\cond
template <typename t1, typename t2>
concept bool compatible_concept = requires (t1, t2)
{
    requires (dimension_v<t1> == dimension_v<t2>);

    requires std::is_same_v<std::decay_t<innermost_value_type_t<t1>>, std::decay_t<innermost_value_type_t<t2>>>;
};
//!\endcond

//!\}

} // namespace seqan3::detail

#ifndef NDEBUG
static_assert(seqan3::has_value_type_v<std::string> == true);
static_assert(seqan3::has_value_type_v<typename std::string::iterator> == true);
static_assert(seqan3::has_value_type_v<typename std::vector<char>::iterator> == true);
static_assert(seqan3::has_value_type_v<int> == false);

static_assert(std::is_same_v<seqan3::innermost_value_type_t<std::string>, char>);
static_assert(std::is_same_v<seqan3::innermost_value_type_t<std::vector<std::string>>, char>);
static_assert(std::is_same_v<seqan3::innermost_value_type_t<std::vector<std::vector<char>>>, char>);

static_assert(seqan3::dimension_v<std::string> == 1);
static_assert(seqan3::dimension_v<std::vector<std::string>> == 2);
static_assert(seqan3::dimension_v<std::vector<std::vector<char>>> == 2);

static_assert(seqan3::compatible_concept<std::string, std::vector<char>>);
static_assert(seqan3::compatible_concept<std::vector<std::string>, std::vector<std::vector<char>>>);
#endif
