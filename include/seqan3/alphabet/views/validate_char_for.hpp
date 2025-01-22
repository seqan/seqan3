// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::views::validate_char_for.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/utility/type_traits/basic.hpp>
#include <seqan3/utility/views/deep.hpp>

namespace seqan3::views
{
/*!\brief               An identity view that throws if an encountered character is not valid for the given alphabet.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \tparam alphabet_t   The alphabet to check validity for; must satisfy seqan3::alphabet.
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             The range of input values. See below for the properties of the returned range.
 * \throws seqan3::invalid_char_assignment if an invalid character is encountered.
 * \ingroup alphabet_views
 *
 * \details
 *
 * \header_file{seqan3/alphabet/views/validate_char_for.hpp}
 *
 * This view throws if it encounters a character that is not in the valid set of an alphabet. It performs
 * no transformation on the elements of the view itself. However, the contiguous property is lost due to the way it
 * is currently implemented.
 *
 * ### View properties
 *
 * This view is a **deep view**. Given a range-of-range as input (as opposed to just a range), it will apply
 * the transformation on the innermost range (instead of the outermost range).
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)         |
 * |----------------------------------|:-------------------------------------:|:--------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                            |
 * | std::ranges::forward_range       |                                       | *preserved*                            |
 * | std::ranges::bidirectional_range |                                       | *preserved*                            |
 * | std::ranges::random_access_range |                                       | *preserved*                            |
 * | std::ranges::contiguous_range    |                                       | *lost*                                 |
 * |                                  |                                       |                                        |
 * | std::ranges::viewable_range      | *required*                            | *guaranteed*                           |
 * | std::ranges::view                |                                       | *guaranteed*                           |
 * | std::ranges::sized_range         |                                       | *preserved*                            |
 * | std::ranges::common_range        |                                       | *preserved*                            |
 * | std::ranges::output_range        |                                       | *preserved*                            |
 * | std::semiregular                 |                                       | *preserved*                            |
 * | seqan3::const_iterable_range     |                                       | *preserved*                            |
 * |                                  |                                       |                                        |
 * | std::ranges::range_reference_t   | seqan3::alphabet_char_t<alphabet_t>   | std::ranges::range_reference_t<urng_t> |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/views/validate_char_for.cpp
 * \hideinitializer
 *
 * \experimentalapi{Experimental since version 3.2.}
 */
template <alphabet alphabet_type>
inline auto const validate_char_for = deep{std::views::transform(
    []<typename char_t>(char_t && in) -> char_t
    {
        static_assert(
            std::common_reference_with<char_t, alphabet_char_t<alphabet_type>>,
            "The innermost value type must have a common reference to underlying char type of alphabet_type.");

        if (!char_is_valid_for<alphabet_type>(in))
        {
            throw seqan3::invalid_char_assignment{"alphabet_type", in};
        }
        return std::forward<char_t>(in);
    })};

} // namespace seqan3::views
