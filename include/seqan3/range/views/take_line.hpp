// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::views::take_line and seqan3::views::take_line_or_throw.
 */

#pragma once

#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/range/views/take_until.hpp>

// ============================================================================
//  views::take_line (adaptor instance definition)
// ============================================================================

namespace seqan3::views
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view adaptor that returns a single line from the underlying range or the full range if there
 *                      is no newline.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             All characters of the underlying range up until, but excluding a unix or windows end-line
 *                      (`\n` or `\r\n`). See below for the properties of the returned range.
 * \ingroup views
 *
 * \details
 *
 * \header_file{seqan3/range/views/take_line.hpp}
 *
 * This adaptor returns a single line **excluding** the end-line character(s), *but moving the cursor behind them
 * for single-pass ranges.* I.e. for all ranges that satisfy std::ranges::forward_range this is the same as calling
 * std::views::take_while:
 * \snippet test/snippet/range/views/take_line_adaptor_def.cpp usage
 * but for *single pass input ranges* this means that any endline characters after the returned range are also consumed
 * (this potentially includes multiple newline characters).
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
 * | std::ranges::sized_range         |                                       | *lost*                                             |
 * | std::ranges::common_range        |                                       | *lost*                                             |
 * | std::ranges::output_range        |                                       | *preserved*                                        |
 * | seqan3::const_iterable_range     |                                       | *preserved*                                        |
 * |                                  |                                       |                                                    |
 * | std::ranges::range_reference_t   | std::common_reference_with<char>      | std::ranges::range_reference_t<urng_t>             |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * Behaviour on std::ranges::forward_range:
 * \include test/snippet/range/views/take_line_behaviour.cpp
 *
 * On single pass std::ranges::input_range it can be used to tokenise the input stream line-wise:
 * \include test/snippet/range/views/take_line_tokenise.cpp
 */
inline auto constexpr take_line = views::take_until_and_consume(is_char<'\r'> || is_char<'\n'>);

// ============================================================================
//  views::take_line_or_throw (adaptor instance definition)
// ============================================================================

/*!\brief A view adaptor that returns a single line from the underlying range (throws if there is no end-of-line).
 * \throws seqan3::unexpected_end_of_input If the underlying range contains no end-of-line marker.
 * \ingroup views
 *
 * \copydetails seqan3::views::take_line
 */
inline auto constexpr take_line_or_throw = views::take_until_or_throw_and_consume(is_char<'\r'> || is_char<'\n'>);

//!\}

} // namespace seqan3::views
