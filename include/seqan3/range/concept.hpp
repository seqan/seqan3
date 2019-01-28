// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Additional non-standard concepts for ranges.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/ranges>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\interface seqan3::const_iterable_concept <>
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
SEQAN3_CONCEPT const_iterable_concept =
    std::ranges::InputRange<std::remove_const_t<type>> &&
    std::ranges::InputRange<type const> &&
    (std::ranges::ForwardRange<std::remove_const_t<type>>       == std::ranges::ForwardRange<type const>) &&
    (std::ranges::BidirectionalRange<std::remove_const_t<type>> == std::ranges::BidirectionalRange<type const>) &&
    (std::ranges::RandomAccessRange<std::remove_const_t<type>>  == std::ranges::RandomAccessRange<type const>);
//!\endcond

} // namespace seqan3
