// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Tobias Loka <LokaT AT rki.de>
 * \brief Provides seqan3::views::to_lower.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/char_operations/transform.hpp>
#include <seqan3/range/views/deep.hpp>
#include <seqan3/std/ranges>

namespace seqan3::views
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view that calls seqan3::to_lower() on each element in the input range.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of converted elements. See below for the properties of the returned range.
 * \ingroup views
 *
 * \details
 *
 * \header_file{seqan3/range/views/to_lower.hpp}
 *
 * ### View properties
 *
 * This view is a **deep view** Given a range-of-range as input (as opposed to just a range), it will apply
 * the transformation on the innermost range (instead of the outermost range).
 *
 * | Concepts and traits              | `urng_t` (underlying range type) | `rrng_t` (returned range type)                                     |
 * |----------------------------------|:--------------------------------:|:------------------------------------------------------------------:|
 * | std::ranges::input_range         | *required*                       | *preserved*                                                        |
 * | std::ranges::forward_range       |                                  | *preserved*                                                        |
 * | std::ranges::bidirectional_range |                                  | *preserved*                                                        |
 * | std::ranges::random_access_range |                                  | *preserved*                                                        |
 * | std::ranges::contiguous_range    |                                  | *lost*                                                             |
 * |                                  |                                  |                                                                    |
 * | std::ranges::viewable_range      | *required*                       | *guaranteed*                                                       |
 * | std::ranges::view                |                                  | *guaranteed*                                                       |
 * | std::ranges::sized_range         |                                  | *preserved*                                                        |
 * | std::ranges::common_range        |                                  | *preserved*                                                        |
 * | std::ranges::output_range        |                                  | *lost*                                                             |
 * | seqan3::const_iterable_range     |                                  | *preserved*                                                        |
 * |                                  |                                  |                                                                    |
 * | std::ranges::range_reference_t   | seqan3::builtin_character        | seqan3::remove_reference_t<std::ranges::range_reference_t<urngt_>> |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 * \include test/snippet/range/views/to_lower.cpp
 * \hideinitializer
 */
inline auto const to_lower = deep{std::views::transform([] (auto const in) noexcept
{
    static_assert(builtin_character<remove_cvref_t<decltype(in)>>,
                  "The value type of seqan3::views::to_lower must model seqan3::builtin_character.");
    return seqan3::to_lower(in);
})};

//!\}

} // namespace seqan3::views
