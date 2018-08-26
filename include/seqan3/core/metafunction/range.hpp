// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
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
#include <seqan3/core/metafunction/pre.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/iterator.hpp>
#include <seqan3/std/ranges>
#include <seqan3/std/iterator>

// TODO(h-2): add innermost_reference instead of or addition to innermost_value_type?

//NOTE(h-2): for the range overloads we explicitly forbid that the type is iteratoer
// because some types are actually both (e.g. std::directory_iterator)

namespace seqan3
{

/*!\addtogroup metafunction
 * \{
 */

// ----------------------------------------------------------------------------
// iterator
// ----------------------------------------------------------------------------

/*!\brief For a seqan3::range return its iterator type. [Type metafunction].
 * \tparam t The type to operate on.
 */
template <std::ranges::Range rng_t>
using iterator_t = ranges::iterator_t<rng_t>;

// ----------------------------------------------------------------------------
// sentinel
// ----------------------------------------------------------------------------

/*!\brief For a seqan3::range return its sentinel type. [Type metafunction].
 * \tparam t The type to operate on.
 */
template <std::ranges::Range rng_t>
using sentinel_t = ranges::sentinel_t<rng_t>;

// ----------------------------------------------------------------------------
// value_type
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `value_type` of another type [specialisation for input ranges].
 * \tparam t The type you wish to query; must satisfy std::ranges::InputRange.
 */
template <std::ranges::InputRange rng_t>
//!\cond
    requires !std::Iterator<rng_t>
//!\endcond
struct value_type<rng_t>
{
    //!\brief Return the value_type member definition from the queried type's iterator.
    using type = value_type_t<iterator_t<rng_t>>;
};

// ----------------------------------------------------------------------------
// reference
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `reference` of another type [specialisation for input ranges].
 * \tparam t The type you wish to query; must satisfy std::ranges::InputRange.
 */
template <std::ranges::InputRange rng_t>
//!\cond
    requires !std::Iterator<rng_t>
//!\endcond
struct reference<rng_t>
{
    //!\brief Return the reference member definition from the queried type's iterator.
    using type = reference_t<iterator_t<rng_t>>;
};

// ----------------------------------------------------------------------------
// rvalue_reference
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `rvalue_reference` of another type [specialisation for input ranges].
 * \tparam t The type you wish to query; must satisfy std::ranges::InputRange.
 */
template <std::ranges::InputRange rng_t>
//!\cond
    requires !std::Iterator<rng_t>
//!\endcond
struct rvalue_reference<rng_t>
{
    //!\brief Return the rvalue_reference member definition from the queried type's iterator.
    using type = rvalue_reference_t<iterator_t<rng_t>>;
};

// ----------------------------------------------------------------------------
// const_reference
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `const_reference` of another type [specialisation for input ranges].
 * \tparam t The type you wish to query; must satisfy std::ranges::InputRange.
 */
template <std::ranges::InputRange rng_t>
//!\cond
    requires !std::Iterator<rng_t>
//!\endcond
struct const_reference<rng_t>
{
    //!\brief Resolves to the reference type of the `const_iterator` of t (not the `const iterator`!).
    using type = reference_t<iterator_t<rng_t const>>;
};

// ----------------------------------------------------------------------------
// difference_type
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `difference_type` of another type [specialisation for ranges].
 * \tparam t The type you wish to query; must satisfy std::ranges::InputRange.
 */
template <std::ranges::Range rng_t>
//!\cond
    requires !std::Iterator<rng_t>
//!\endcond
struct difference_type<rng_t>
{
    //!\brief Return the difference_type member definition from the queried type's iterator.
    using type = difference_type_t<iterator_t<rng_t>>;
};

// ----------------------------------------------------------------------------
// size_type
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `size_type` of another type [specialisation for sized ranges].
 * \tparam t The type you wish to query; must satisfy std::ranges::SizedRange.
 */
template <std::ranges::SizedRange rng_t>
//!\cond
    requires !std::Iterator<rng_t>
//!\endcond
struct size_type<rng_t>
{
    //!\brief Return the size_type member definition from the queried type's iterator.
    using type = size_type_t<iterator_t<rng_t>>;
};

// ----------------------------------------------------------------------------
// innermost_value_type
// ----------------------------------------------------------------------------

//NOTE(h-2): this could be moved to a separate file, because it also applies to iterators

/*!\brief Recursively determines the `value_type` on containers and/or iterators [Type metafunction].
 * \tparam t The type to recurse on; must have `ranges::value_type_t<rng_t>`
 *
 * \details
 *
 * Attention, this metafunction implicitly removes cv-qualifiers on the all value_types except the one returned.
 */
template <typename t>
//!\cond
    requires requires (t) { typename value_type_t<remove_cvref_t<t>>; }
//!\endcond
struct innermost_value_type
{
    //!\brief The return type (recursion not shown).
    using type = value_type_t<remove_cvref_t<t>>;
};

//!\cond
template <typename t>
    requires requires (t)
    {
        typename value_type_t<remove_cvref_t<t>>;
        typename value_type_t<remove_cvref_t<value_type_t<remove_cvref_t<t>>>>;
    }
struct innermost_value_type<t>
{
    using type = typename innermost_value_type<value_type_t<remove_cvref_t<t>>>::type;
};
//!\endcond

//!\brief Shortcut for seqan3::innermost_value_type.
template <typename t>
using innermost_value_type_t = typename innermost_value_type<t>::type;

// ----------------------------------------------------------------------------
// dimension_v
// ----------------------------------------------------------------------------

//NOTE(h-2): this could be moved to a separate file, because it also applies to iterators

/*!\brief Returns the number of times you can call `seqan3::value_type_t` recursively on t [Value metafunction].
 * \tparam t The type to be queried; must resolve `seqan3::value_type_t` at least once.
 *
 * \details
 *
 * Attention, this metafunction implicitly removes cv-qualifiers and reference from the types it recurses on and
 * returns.
 */
template <typename t>
//!\cond
    requires requires (t) { typename value_type_t<remove_cvref_t<t>>; }
//!\endcond
constexpr size_t dimension_v = 1;

//!\cond
template <typename t>
    requires requires (t)
    {
        typename value_type_t<remove_cvref_t<t>>;
        typename value_type_t<value_type_t<remove_cvref_t<t>>>;
    }
constexpr size_t dimension_v<t> = dimension_v<value_type_t<remove_cvref_t<t>>> + 1;
//!\endcond

// ----------------------------------------------------------------------------
// compatible_concept
// ----------------------------------------------------------------------------

//NOTE(h-2): this could be moved to a separate file, because it also applies to iterators

/*!\interface seqan3::compatible_concept <>
 * \brief Two types are "compatible" if their seqan3::dimension_v and their seqan3::innermost_value_type_t are
 * the same.
 *
 * \details
 *
 * ```cpp
 * // these evaluate to true:
 * static_assert(seqan3::compatible_concept<std::string,              std::vector<char>>);
 * static_assert(seqan3::compatible_concept<std::vector<std::string>, std::vector<std::vector<char>>>);
 * ```
 *
 * Attention, this metafunction implicitly removes cv-qualifiers and reference from the types it recurses on and
 * compares.
 */
//!\cond
template <typename t1, typename t2>
concept bool compatible_concept = requires (t1, t2)
{
    requires (dimension_v<t1> == dimension_v<t2>);

    requires std::is_same_v<innermost_value_type_t<t1>, innermost_value_type_t<t2>>;
};
//!\endcond

//!\}

} // namespace seqan3
