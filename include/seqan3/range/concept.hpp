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
 * \brief Adaptations of concepts from the Ranges TS
 * \ingroup range
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <range/v3/range_concepts.hpp>

namespace seqan3
{
/*!\addtogroup range
 * \{
 */
/*!\interface seqan3::range_concept <>
 * \brief Defines the requirements of a type that allows iteration over its elements by providing a begin iterator
 * and an end sentinel.
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/Range
 */
//!\cond
template <typename type>
concept bool range_concept               = (bool)ranges::Range<type>();
//!\endcond

/*!\interface seqan3::sized_range_concept <>
 * \extends seqan3::range_concept
 * \brief Specifies the requirements of a Range type that knows its size in constant time with the size function.
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/SizedRange
 */
//!\cond
template <typename type>
concept bool sized_range_concept         = range_concept<type> && (bool)ranges::SizedRange<type>();
//!\endcond

/*!\interface seqan3::bounded_range_concept <>
 * \extends seqan3::range_concept
 * \brief Specifies requirements of a Range type for which `begin` and `end` return objects of the same type.
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/BoundedRange
 */
//!\cond
template <typename type>
concept bool bounded_range_concept       = range_concept<type> && (bool)ranges::BoundedRange<type>();
//!\endcond

/*!\interface seqan3::output_range_concept <>
 * \extends seqan3::range_concept
 * \brief Specifies requirements of a Range type for which `begin` returns a type that satisfies
 * seqan3::output_iterator_concept.
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/OutputRange
 */
//!\cond
template <typename type, typename out_type>
concept bool output_range_concept        = range_concept<type> && (bool)ranges::OutputRange<type, out_type>();
//!\endcond

/*!\interface seqan3::input_range_concept <>
 * \extends seqan3::range_concept
 * \brief Specifies requirements of an Range type for which `begin` returns a type that satisfies
 * seqan3::input_iterator_concept.
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/InputRange
 */
//!\cond
template <typename type>
concept bool input_range_concept         = range_concept<type> && (bool)ranges::InputRange<type>();
//!\endcond

/*!\interface seqan3::forward_range_concept <>
 * \extends seqan3::input_range_concept
 * \brief Specifies requirements of an Range type for which `begin` returns a type that satisfies
 * seqan3::forward_iterator_concept.
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/ForwardRange
 */
//!\cond
template <typename type>
concept bool forward_range_concept       = input_range_concept<type> && (bool)ranges::ForwardRange<type>();
//!\endcond

/*!\interface seqan3::bidirectional_range_concept <>
 * \extends seqan3::forward_range_concept
 * \brief Specifies requirements of an Range type for which `begin` returns a type that satisfies
 * seqan3::bidirectional_iterator_concept.
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/BidirectionalRange
 */
//!\cond
template <typename type>
concept bool bidirectional_range_concept = forward_range_concept<type> && (bool)ranges::BidirectionalRange<type>();
//!\endcond

/*!\interface seqan3::random_access_range_concept <>
 * \extends seqan3::bidirectional_range_concept
 * \brief Specifies requirements of an Range type for which `begin` returns a type that satisfies
 * seqan3::random_access_iterator_concept.
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/RandomAccessRange
 */
//!\cond
template <typename type>
concept bool random_access_range_concept = bidirectional_range_concept<type> && (bool)ranges::RandomAccessRange<type>();
//!\endcond

//!\}
} // namespace seqan3

#ifndef NDEBUG
/* Check the STL containers */

#include <vector>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <string>

// no fwd list
static_assert(seqan3::sized_range_concept<std::list<char>>);
static_assert(seqan3::sized_range_concept<std::array<char, 2>>);
static_assert(seqan3::sized_range_concept<std::vector<char>>);
static_assert(seqan3::sized_range_concept<std::deque<char>>);
static_assert(seqan3::sized_range_concept<std::string>);

static_assert(seqan3::bounded_range_concept<std::forward_list<char>>);
static_assert(seqan3::bounded_range_concept<std::list<char>>);
static_assert(seqan3::bounded_range_concept<std::array<char, 2>>);
static_assert(seqan3::bounded_range_concept<std::vector<char>>);
static_assert(seqan3::bounded_range_concept<std::deque<char>>);
static_assert(seqan3::bounded_range_concept<std::string>);

static_assert(seqan3::output_range_concept<std::forward_list<char>, char>);
static_assert(seqan3::output_range_concept<std::list<char>, char>);
static_assert(seqan3::output_range_concept<std::array<char, 2>, char>);
static_assert(seqan3::output_range_concept<std::vector<char>, char>);
static_assert(seqan3::output_range_concept<std::deque<char>, char>);
static_assert(seqan3::output_range_concept<std::string, char>);

static_assert(seqan3::forward_range_concept<std::forward_list<char>>);
static_assert(seqan3::bidirectional_range_concept<std::list<char>>);
static_assert(seqan3::random_access_range_concept<std::array<char, 2>>);
static_assert(seqan3::random_access_range_concept<std::vector<char>>);
static_assert(seqan3::random_access_range_concept<std::deque<char>>);
static_assert(seqan3::random_access_range_concept<std::string>);

/* Check the SDSL containers */
#include <sdsl/int_vector.hpp>
static_assert(seqan3::sized_range_concept<sdsl::int_vector<>>);
static_assert(seqan3::bounded_range_concept<sdsl::int_vector<>>);
// doesn't work?
// static_assert(seqan3::output_range_concept<sdsl::int_vector<>, int>);
static_assert(seqan3::random_access_range_concept<sdsl::int_vector<>>);

/* Check range-v3 containers */
#include <range/v3/view/any_view.hpp>

static_assert(seqan3::range_concept<ranges::any_random_access_view<char>>);
// static_assert(seqan3::sized_range_concept<ranges::any_random_access_view<char>>);
static_assert(seqan3::random_access_range_concept<ranges::any_random_access_view<char>>);

/* Check our containers */
//TODO
#endif
