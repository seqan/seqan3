// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::views::to_char.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/utility/views/deep.hpp>

namespace seqan3::views
{
/*!\brief               A view that calls seqan3::to_char() on each element in the input range.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of converted elements. See below for the properties of the returned range.
 * \ingroup alphabet_views
 *
 * \details
 *
 * \header_file{seqan3/alphabet/views/to_char.hpp}
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
 * | std::ranges::range_reference_t   | seqan3::alphabet                      | seqan3::alphabet_char_t<std::ranges::range_value_t<urng_t>>   |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 * \include test/snippet/alphabet/views/range_view_to_char.cpp
 * \hideinitializer
 *
 * \stableapi{Since version 3.1.}
 */
inline auto const to_char = deep{std::views::transform(
    [](auto const in) noexcept
    {
        static_assert(alphabet<std::remove_cvref_t<decltype(in)>>,
                      "The value type of seqan3::views::to_char must model the seqan3::alphabet.");
        return seqan3::to_char(in);
    })};

} // namespace seqan3::views
