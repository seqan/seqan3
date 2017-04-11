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

#pragma once

#include <range/v3/range_concepts.hpp>

/*!\file core/concept/range.hpp
 * \brief Adaptations of concepts from the Ranges TS
 * \ingroup core
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

namespace seqan3
{

/*! resolves to `ranges::Range<type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/Range
 */
template <typename type>
concept bool range_concept               = (bool)ranges::Range<type>();

/* these are independent specializations of range_concept */

/*! resolves to `ranges::View<type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/View
 */
template <typename type>
concept bool view_concept                = range_concept<type> && (bool)ranges::View<type>();

/*! resolves to `ranges::SizedRange<type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/SizedRange
 */
template <typename type>
concept bool sized_range_concept         = range_concept<type> && (bool)ranges::SizedRange<type>();

/*! resolves to `ranges::BoundedRange<type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/BoundedRange
 */
template <typename type>
concept bool bounded_range_concept       = range_concept<type> && (bool)ranges::BoundedRange<type>();

/*! resolves to `ranges::OutputRange<type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/OutputRange
 */
template <typename type, typename out_type>
concept bool output_range_concept        = range_concept<type> && (bool)ranges::OutputRange<type, out_type>();


/* the following specialize incrementally */

/*! resolves to `ranges::InputRange<type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/InputRange
 */
template <typename type>
concept bool input_range_concept         = range_concept<type> && (bool)ranges::InputRange<type>();

/*! resolves to `ranges::ForwardRange<type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/ForwardRange
 */
template <typename type>
concept bool forward_range_concept       = input_range_concept<type> && (bool)ranges::ForwardRange<type>();

/*! resolves to `ranges::BidirectionalRange<type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/BidirectionalRange
 */
template <typename type>
concept bool bidirectional_range_concept = forward_range_concept<type> && (bool)ranges::BidirectionalRange<type>();

/*! resolves to `ranges::RandomAccessRange<type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/RandomAccessRange
 */
template <typename type>
concept bool random_access_range_concept = bidirectional_range_concept<type> && (bool)ranges::RandomAccessRange<type>();

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
