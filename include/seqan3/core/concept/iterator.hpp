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

/*!\file core/concept/iterator.hpp
 * \brief Adaptions of Iterator concepts from the Ranges TS.
 * \ingroup core
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <range/v3/range_concepts.hpp>
#include <range/v3/utility/iterator.hpp>

#include <seqan3/core/concepts/core.hpp>

namespace seqan3
{

//!\name Iterator Concepts
//!\{

/*!\brief Resolves to `ranges::Readable<type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/Readable
 */
template <typename t>
concept bool readable_concept =                 static_cast<bool>(ranges::Readable<t>());

/*!\brief Resolves to `ranges::Writable<out_type, type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/Writable
 */
template <typename out, typename t>
concept bool writable_concept =                 static_cast<bool>(ranges::Writable<out, t>());

/*!\brief Resolves to `ranges::WeaklyIncrementable<type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/WeaklyIncrementable
 */
template<typename i>
concept bool weakly_incrementable_concept =     semi_regular_concept<i> &&
                                                static_cast<bool>(ranges::WeaklyIncrementable<i>());
/*!\brief Resolves to `ranges::Incrementable<type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/Incrementable
 */
template<typename i>
concept bool incrementable_concept =            regular_concept<i> &&
                                                weakly_incrementable_concept<i> &&
                                                static_cast<bool>(ranges::Incrementable<i>());
/*!\brief Resolves to `ranges::Iterator<iterator_type>()`
 * \sa http://en.cppreference.com/w/cpp/concept/Iterator
 */
template<typename i>
concept bool iterator_concept =                 weakly_incrementable_concept<i> &&
                                                copyable_concept<i> &&
                                                static_cast<bool>(ranges::Iterator<i>());

/*!\brief Resolves to `ranges::Sentinel<sentinel_type, iterator_type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/Sentinel
 */
template<typename s, typename i>
concept bool sentinel_concept =                 semi_regular_concept<s> &&
                                                iterator_concept<i> &&
                                                static_cast<bool>(ranges::Sentinel<s, i>());

/*!\brief Resolves to `ranges::SizedSentinel<sentinel_type, iterator_type>()`
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/SizedSentinel
 */
template<typename s, typename i>
concept bool sized_sentinel_concept =           sentinel_concept<s, i> &&
                                                static_cast<bool>(ranges::SizedSentinel<s, i>());

/*!\brief Resolves to `ranges::OutputIterator<iterator_type, type>()`
 * \sa http://en.cppreference.com/w/cpp/concept/OutputIterator
 */
template<typename out, typename t>
concept bool output_iterator_concept =          iterator_concept<out> &&
                                                writable_concept<out, t> &&
                                                static_cast<bool>(ranges::OutputIterator<out, t>());

/*!\brief Resolves to `ranges::InputIterator<iterator_type>()`
 * \sa http://en.cppreference.com/w/cpp/concept/InputIterator
 */
template<typename i>
concept bool input_iterator_concept =           iterator_concept<i> &&
                                                readable_concept<i> &&
                                                static_cast<bool>(ranges::InputIterator<i>());

/*!\brief Resolves to `ranges::ForwardIterator<iterator_type>()`
 * \sa http://en.cppreference.com/w/cpp/concept/ForwardIterator
 */
template<typename i>
concept bool forward_iterator_concept =         input_iterator_concept<i> &&
                                                incrementable_concept<i> &&
                                                sentinel_concept<i, i> &&
                                                static_cast<bool>(ranges::ForwardIterator<i>());

/*!\brief Resolves to `ranges::BidirectionalIterator<iterator_type>()`
 * \sa http://en.cppreference.com/w/cpp/concept/BidirectionalIterator
 */
template<typename i>
concept bool bidirectional_iterator_concept =   forward_iterator_concept<i> &&
                                                static_cast<bool>(ranges::BidirectionalIterator<i>());

/*!\brief Resolves to `ranges::RandomAccessIterator<iterator_type>()`
 * \sa http://en.cppreference.com/w/cpp/concept/RandomAccessIterator
 */
template<typename i>
concept bool random_access_iterator_concept =   bidirectional_iterator_concept<i> &&
                                                totally_ordered_concept<i> &&
                                                sized_sentinel_concept<i, i> &&
                                                static_cast<bool>(ranges::RandomAccessIterator<i>());
//!\} // Iterator Concepts.
}  // namespace seqan3

#ifndef NDEBUG
/* Check the iterator concepts */

#include <seqan3/core/concepts/iterator_detail.hpp>

#endif // NDEBUG
