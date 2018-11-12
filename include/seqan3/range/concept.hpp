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
 * \brief Additional non-standard concepts for ranges.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/ranges>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\interface seqan3::ConstIterableRange <>
 * \extends std::InputRange
 * \brief Specifies requirements of an input range type for which the `const` version of that type satisfies the
 * same strength input range concept as the non-const version.
 *
 * \details
 *
 * For a type `t` it usually holds that if `t` is a range, `t const` is also a range with similar properties, but
 * there are cases where this does not hold:
 *
 *   * a `const` range is usually not writable so std::OutputRange is lost; pure output ranges
 * (those that are not also input ranges) are therefore not `const`-iterable;
 *   * single-pass input ranges, like SeqAn files, are not `const`-iterable, because "single-pass-ness" implies that
 * there is something in the range that changes on every iterator increment (and `const` ranges can't change);
 *   * certain views store a state with their algorithm that also changes when `begin()` is called or an
 * iterator is incremented; these may be not be `const`-iterable, because the standard library
 * (and also SeqAn3) guarantees that it is safe to call `const`-qualified functions concurrently.
 */
//!\cond
template <typename type>
concept ConstIterableRange =
    std::ranges::InputRange<std::remove_const_t<type>> &&
    std::ranges::InputRange<type const> &&
    (std::ranges::ForwardRange<std::remove_const_t<type>>       == std::ranges::ForwardRange<type const>) &&
    (std::ranges::BidirectionalRange<std::remove_const_t<type>> == std::ranges::BidirectionalRange<type const>) &&
    (std::ranges::RandomAccessRange<std::remove_const_t<type>>  == std::ranges::RandomAccessRange<type const>);
//!\endcond

} // namespace seqan3
