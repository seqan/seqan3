// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::views::take_exactly and seqan3::views::take_exactly_or_throw.
 */

#pragma once

#include <seqan3/range/views/take.hpp>

// ============================================================================
//  views::take_exactly (adaptor instance definition)
// ============================================================================

namespace seqan3::views
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
 * \ingroup views
 *
 * \details
 *
 * \header_file{seqan3/range/views/take_exactly.hpp}
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |----------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                                        |
 * | std::ranges::forward_range       |                                       | *preserved*                                        |
 * | std::ranges::bidirectional_range |                                       | *preserved*                                        |
 * | std::ranges::random_access_range |                                       | *preserved*                                        |
 * | std::ranges::contiguous_range    |                                       | *preserved*                                        |
 * |                                  |                                       |                                                    |
 * | std::ranges::viewable_range      | *required*                            | *guaranteed*                                       |
 * | std::ranges::view                |                                       | *guaranteed*                                       |
 * | std::ranges::sized_range         |                                       | ***guaranteed***                                   |
 * | std::ranges::common_range        |                                       | *preserved*                                        |
 * | std::ranges::output_range        |                                       | *preserved* except if `urng_t` is std::basic_string|
 * | seqan3::const_iterable_range     |                                       | *preserved*                                        |
 * |                                  |                                       |                                                    |
 * | std::ranges::range_reference_t   |                                       | std::ranges::range_reference_t<urng_t>             |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * The difference to seqan3::views::take is that this view always exposes size information – even if the underlying
 * range is not sized. You should only use this if you know that the underlying range will always be
 * at least `size` long.
 *
 * For seqan3::views::take_exactly if the underlying range is shorter than `size`, the behaviour is undefined.
 * seqan3::views::take_exactly_or_throw is a safer alternative, because it throws an exception when an iterator before
 * the `size`-th one compares equal to the end sentinel; and it also throws on construction if it knows that the
 * underlying range is smaller.
 *
 * ### Example
 *
 * \include test/snippet/range/views/take_exactly.cpp
 *
 * \hideinitializer
 */
inline auto constexpr take_exactly = detail::take_fn<true, false>{};

// ============================================================================
//  views::take_exactly_or_throw (adaptor instance definition)
// ============================================================================

/*!\brief A view adaptor that returns the first `size` elements from the underlying range and also exposes size
 *        information; throws if the underlying range is smaller than `size`.
 * \throws seqan3::unexpected_end_of_input If the underlying range is smaller than `size`.
 * \ingroup views
 *
 * \copydetails seqan3::views::take_exactly
 * \hideinitializer
 */
inline auto constexpr take_exactly_or_throw = detail::take_fn<true, true>{};

//!\}
} // namespace seqan3::views
