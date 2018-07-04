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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <range/v3/range_concepts.hpp>
#include <seqan3/core/platform.hpp>

namespace seqan3
{
/*!\addtogroup concept
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

/*!\interface seqan3::common_range_concept  <>
 * \extends seqan3::range_concept
 * \brief Specifies requirements of a Range type for which `begin` and `end` return objects of the same type.
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/BoundedRange
 */
//!\cond
template <typename type>
concept bool common_range_concept        = range_concept<type> && (bool)ranges::BoundedRange<type>();
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

/*!\interface seqan3::contiguous_range_concept <>
 * \extends seqan3::random_access_range_concept
 * \brief Specifies requirements of an Range type whose elements occupy adjacent locations in memory.
 */
//!\cond
template <typename type>
concept bool contiguous_range_concept = random_access_range_concept<type> && (bool)ranges::ContiguousRange<type>();
//!\endcond

/*!\interface seqan3::const_iterable_concept <>
 * \extends seqan3::input_range_concept
 * \brief Specifies requirements of an input range type for which the `const` version of that type satisfies the
 * same strength input range concept as the non-const version.
 *
 * \details
 *
 * For a type `t` it usually holds that if `t` is a range, `t const` is also a range with similar properties, but
 * there are cases where this does not hold:
 *
 *   * a `const` range is usually not writable so seqan3::output_range_concept is lost; pure output ranges
 * (those that are not also input ranges) are therefore not `const`-iterable;
 *   * single-pass input ranges, like SeqAn files, are not `const`-iterable, because "single-pass-ness" implies that
 * there is something in the range that changes on every iterator increment (and `const` ranges can't change);
 *   * certain views store a state with their algorithm that also changes when `begin()` is called or an
 * iterator is incremented; these may be not be `const`-iterable, because the standard library
 * (and also SeqAn3) guarantees that it is safe to call `const`-qualified functions concurrently.
 */
//!\cond
template <typename type>
concept bool const_iterable_concept =
    input_range_concept<std::remove_const_t<type>> &&
    input_range_concept<type const> &&
    (forward_range_concept<std::remove_const_t<type>>       == forward_range_concept<type const>) &&
    (bidirectional_range_concept<std::remove_const_t<type>> == bidirectional_range_concept<type const>) &&
    (random_access_range_concept<std::remove_const_t<type>> == random_access_range_concept<type const>);
//!\endcond


/*!\interface seqan3::view_concept <>
 * \extends seqan3::semi_regular_concept
 * \extends seqan3::viewable_range_concept
 * \brief Specifies the requirements of a Range type that has constant time copy, move and assignment operators.
 * \sa \ref view
 * \sa http://en.cppreference.com/w/cpp/experimental/ranges/iterator/View
 */
//!\cond
template <typename type>
concept bool view_concept = range_concept<type> && (bool)ranges::View<type>();
//!\endcond

/*!\interface seqan3::viewable_range_concept <>
 * \extends seqan3::range_concept
 * \brief Specifies the requirements of a Range type that is either a seqan3::view_concept or an lvalue-reference.
 * \sa \ref view
 */
//!\cond
template <typename type>
concept bool viewable_range_concept = range_concept<type> && (std::is_lvalue_reference_v<type> || view_concept<type>);
//!\endcond

//!\}
} // namespace seqan3
