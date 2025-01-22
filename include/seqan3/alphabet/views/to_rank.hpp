// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::views::to_rank.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cstring>
#include <ranges>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/utility/views/deep.hpp>

namespace seqan3::views
{
/*!\brief               A view that calls seqan3::to_rank() on each element in the input range.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of converted elements. See below for the properties of the returned range.
 * \ingroup alphabet_views
 *
 * \details
 *
 * \header_file{seqan3/alphabet/views/to_rank.hpp}
 *
 * ### View properties
 *
 * This view is a **deep view** Given a range-of-range as input (as opposed to just a range), it will apply
 * the transformation on the innermost range (instead of the outermost range).
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                                |
 * |----------------------------------|:-------------------------------------:|:-------------------------------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                                                   |
 * | std::ranges::forward_range       |                                       | *preserved*                                                   |
 * | std::ranges::bidirectional_range |                                       | *preserved*                                                   |
 * | std::ranges::random_access_range |                                       | *preserved*                                                   |
 * | std::ranges::contiguous_range    |                                       | *lost*                                                        |
 * |                                  |                                       |                                                               |
 * | std::ranges::viewable_range      | *required*                            | *guaranteed*                                                  |
 * | std::ranges::view                |                                       | *guaranteed*                                                  |
 * | std::ranges::sized_range         |                                       | *preserved*                                                   |
 * | std::ranges::common_range        |                                       | *preserved*                                                   |
 * | std::ranges::output_range        |                                       | *lost*                                                        |
 * | seqan3::const_iterable_range     |                                       | *preserved*                                                   |
 * |                                  |                                       |                                                               |
 * | std::ranges::range_reference_t   | seqan3::alphabet                      | seqan3::alphabet_rank_t<std::ranges::range_value_t<urng_t>>   |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 * \include test/snippet/alphabet/views/range_view_to_rank.cpp
 * We also convert to unsigned here, because the seqan3::alphabet_rank_t is often `uint8_t` which is
 * often implemented as `unsigned char` and thus will not be printed as a number by default.
 * \hideinitializer
 *
 * \stableapi{Since version 3.1.}
 */
inline auto const to_rank = deep{std::views::transform(
    [](auto const in) noexcept
    {
        static_assert(semialphabet<decltype(in)>,
                      "The value type of seqan3::views::to_rank must model the seqan3::alphabet.");
        return seqan3::to_rank(in);
    })};

} // namespace seqan3::views
