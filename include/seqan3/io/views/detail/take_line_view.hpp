// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::detail::take_line and seqan3::detail::take_line_or_throw.
 */

#pragma once

#include <seqan3/io/views/detail/take_until_view.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

// ============================================================================
//  detail::take_line (adaptor instance definition)
// ============================================================================

namespace seqan3::detail
{
/*!\brief               A view adaptor that returns a single line from the underlying range or the full range if there
 *                      is no newline.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             All characters of the underlying range up until, but excluding a unix or windows end-line
 *                      (`\n` or `\r\n`). See below for the properties of the returned range.
 * \ingroup io_views
 *
 * \details
 *
 * \header_file{seqan3/io/views/detail/take_line_view.hpp}
 *
 * This adaptor returns a single line **excluding** the end-line character(s), *but moving the cursor behind them
 * for single-pass ranges.* I.e. for all ranges that satisfy std::ranges::forward_range this is the same as calling
 * std::views::take_while:
 * \snippet test/snippet/io/views/detail/take_line_view_adaptor_def.cpp usage
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
 * \include test/snippet/io/views/detail/take_line_view_behaviour.cpp
 *
 * On single pass std::ranges::input_range it can be used to tokenise the input stream line-wise:
 * \include test/snippet/io/views/detail/take_line_view_tokenise.cpp
 */
inline constexpr auto take_line = detail::take_until_and_consume(is_char<'\r'> || is_char<'\n'>);

// ============================================================================
//  detail::take_line_or_throw (adaptor instance definition)
// ============================================================================

/*!\brief A view adaptor that returns a single line from the underlying range (throws if there is no end-of-line).
 * \throws seqan3::unexpected_end_of_input If the underlying range contains no end-of-line marker.
 * \ingroup io_views
 *
 * \copydetails seqan3::detail::take_line
 */
inline constexpr auto take_line_or_throw = detail::take_until_or_throw_and_consume(is_char<'\r'> || is_char<'\n'>);

} // namespace seqan3::detail
