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
 * \brief Provides various metafunctions for use on iterators.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <type_traits>

#include <range/v3/utility/iterator_traits.hpp>

#include <seqan3/core/platform.hpp>
#include <seqan3/core/metafunction/pre.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

namespace seqan3
{

/*!\addtogroup metafunction
 * \{
 */

// ----------------------------------------------------------------------------
// value_type
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `value_type` of another type [specialisation for input iterators].
 * \tparam it_t The type you wish to query; must satisfy std::input_Iterator.
 */
template <std::InputIterator it_t>
struct value_type<it_t>
{
    //!\brief Return the member type as return type.
    using type = typename std::iterator_traits<std::remove_reference_t<it_t>>::value_type;
};

// see specialisation for ranges in core/metafunction/range.hpp

// ----------------------------------------------------------------------------
// reference
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `reference` of another type [specialisation for input iterators].
 * \tparam it_t The type you wish to query; must satisfy std::input_Iterator.
 */
template <std::InputIterator it_t>
struct reference<it_t>
{
    //!\brief Return the member type as return type.
    using type = typename std::iterator_traits<std::remove_reference_t<it_t>>::reference;
};

// see specialisation for ranges in core/metafunction/range.hpp

// ----------------------------------------------------------------------------
// rvalue_reference
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `rvalue_reference` of another type [specialisation for input iterators].
 * \tparam it_t The type you wish to query; must satisfy std::input_Iterator.
 */
template <std::InputIterator it_t>
struct rvalue_reference<it_t>
{
    //!\brief Return the member type as return type.
    using type = decltype(std::ranges::iter_move(std::declval<it_t &>()));
};

// see specialisation for ranges in core/metafunction/range.hpp

// ----------------------------------------------------------------------------
// const_reference
// ----------------------------------------------------------------------------

// only defined for ranges

// ----------------------------------------------------------------------------
// difference_type
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `difference_type` of another type [specialisation for iterators].
 * \tparam it_t The type you wish to query; must satisfy std::WeaklyIncrementable.
 */
template <std::WeaklyIncrementable it_t>
struct difference_type<it_t>
{
    //!\brief Return the member type as return type.
    using type = typename std::iterator_traits<std::remove_reference_t<it_t>>::difference_type;
};

// see specialisation for ranges in core/metafunction/range.hpp

// ----------------------------------------------------------------------------
// size_type
// ----------------------------------------------------------------------------

/*!\brief Type metafunction that returns the `size_type` of another type [specialisation for iterators].
 * \tparam it_t The type you wish to query; must satisfy std::WeaklyIncrementable.
 */
template <std::WeaklyIncrementable it_t>
struct size_type<it_t>
{
    //!\brief Return the member type as return type.
    using type = std::make_unsigned_t<difference_type_t<it_t>>;
};

// see specialisation for ranges in core/metafunction/range.hpp

//!\}

} // namespace seqan3
