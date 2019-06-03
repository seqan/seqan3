// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::take_exactly and seqan3::view::take_exactly_or_throw.
 */

#pragma once

#include <seqan3/range/view/take.hpp>

// ============================================================================
//  view::take_exactly (adaptor instance definition)
// ============================================================================

namespace seqan3::view
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view adaptor that returns the first `size` elements from the underlying range (or less if the
 *                      underlying range is shorter); also provides size information.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \param[in] size      The target size of the view.
 * \returns             Up to `size` elements of the underlying range.
 * \ingroup view
 *
 * \details
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/view/take_exactly.hpp>
 * ```
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |---------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                        |
 * | std::ranges::ForwardRange       |                                       | *preserved*                                        |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                                        |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                                        |
 * | std::ranges::ContiguousRange    |                                       | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                       |
 * | std::ranges::View               |                                       | *guaranteed*                                       |
 * | std::ranges::SizedRange         |                                       | ***guaranteed***                                   |
 * | std::ranges::CommonRange        |                                       | *preserved*                                        |
 * | std::ranges::OutputRange        |                                       | *preserved* except if `urng_t` is std::basic_string|
 * | seqan3::ConstIterableRange      |                                       | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | seqan3::reference_t             |                                       | seqan3::reference_t<urng_t>                        |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * The difference to seqan3::view::take is that this view always exposes size information – even if the underlying
 * range is not sized. You should only use this if you know that the underlying range will always be
 * at least `size` long.
 *
 * For seqan3::view::take_exactly if the underlying range is shorter than `size`, the behaviour is undefined.
 * seqan3::view::take_exactly_or_throw is a safer alternative, because it throws an exception when an iterator before
 * the `size`-th one compares equal to the end sentinel; and it also throws on construction if it knows that the
 * underlying range is smaller.
 *
 * ### Example
 *
 * \include test/snippet/range/view/take_exactly.cpp
 *
 * \hideinitializer
 */
inline auto constexpr take_exactly = detail::take_fn<true, false>{};

// ============================================================================
//  view::take_exactly_or_throw (adaptor instance definition)
// ============================================================================

/*!\brief A view adaptor that returns the first `size` elements from the underlying range and also exposes size
 *        information; throws if the underlying range is smaller than `size`.
 * \throws seqan3::unexpected_end_of_input If the underlying range is smaller than `size`.
 * \ingroup view
 *
 * \copydetails seqan3::view::take_exactly
 * \hideinitializer
 */
inline auto constexpr take_exactly_or_throw = detail::take_fn<true, true>{};

//!\}
} // namespace seqan3::view
