// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::views::char_strictly_to.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/validate_char_for.hpp>

namespace seqan3::views
{
/*!\brief               A view over an alphabet, given a range of characters.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \tparam alphabet_t   The alphabet to convert to; must satisfy seqan3::alphabet.
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of converted elements. See below for the properties of the returned range.
 * \throws seqan3::invalid_char_assignment if an invalid character is encountered.
 * \ingroup alphabet_views
 *
 * \details
 *
 * \header_file{seqan3/alphabet/views/char_strictly_to.hpp}
 *
 * This view differs from seqan3::views::chars_to in that it throws an exception if an invalid character conversion
 * happens. See seqan3::char_strictly_to for more details.
 *
 * ### View properties
 *
 * This view is a **deep view**. Given a range-of-range as input (as opposed to just a range), it will apply
 * the transformation on the innermost range (instead of the outermost range).
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type) |
 * |----------------------------------|:-------------------------------------:|:------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                    |
 * | std::ranges::forward_range       |                                       | *preserved*                    |
 * | std::ranges::bidirectional_range |                                       | *preserved*                    |
 * | std::ranges::random_access_range |                                       | *preserved*                    |
 * | std::ranges::contiguous_range    |                                       | *lost*                         |
 * |                                  |                                       |                                |
 * | std::ranges::viewable_range      | *required*                            | *guaranteed*                   |
 * | std::ranges::view                |                                       | *guaranteed*                   |
 * | std::ranges::sized_range         |                                       | *preserved*                    |
 * | std::ranges::common_range        |                                       | *preserved*                    |
 * | std::ranges::output_range        |                                       | *lost*                         |
 * | std::semiregular                 |                                       | *preserved*                    |
 * | seqan3::const_iterable_range     |                                       | *preserved*                    |
 * |                                  |                                       |                                |
 * | std::ranges::range_reference_t   | seqan3::alphabet_char_t<alphabet_t>   | `alphabet_t`                   |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/views/char_strictly_to.cpp
 * \hideinitializer
 *
 * \experimentalapi{Experimental since version 3.2.}
 */
template <alphabet alphabet_type>
inline auto const char_strictly_to = views::validate_char_for<alphabet_type> | views::char_to<alphabet_type>;

} // namespace seqan3::views
