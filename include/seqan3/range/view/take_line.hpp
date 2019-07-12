// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::take_line and seqan3::view::take_line_or_throw.
 */

#pragma once

#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/range/view/take_until.hpp>

// ============================================================================
//  view::take_line (adaptor instance definition)
// ============================================================================

namespace seqan3::view
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
 * \ingroup view
 *
 * \details
 *
 * This adaptor returns a single line **excluding** the end-line character(s), *but moving the cursor behind them
 * for single-pass ranges.* I.e. for all ranges that satisfy std::ranges::ForwardRange this is the same as calling
 * \snippet test/snippet/range/view/take_line.cpp adaptor_def
 * but for *single pass input ranges* this means that any endline characters after the returned range are also consumed
 * (this potentially includes multiple newline characters).
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/view/take_line.hpp>
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
 * | std::ranges::SizedRange         |                                       | *lost*                                             |
 * | std::ranges::CommonRange        |                                       | *lost*                                             |
 * | std::ranges::OutputRange        |                                       | *preserved*                                        |
 * | seqan3::ConstIterableRange      |                                       | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | seqan3::reference_t             | std::CommonReference<char>            | seqan3::reference_t<urng_t>                        |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * Behaviour on std::ranges::ForwardRange:
 * \snippet test/snippet/range/view/take_line.cpp behaviour
 *
 * On single pass std::ranges::InputRange it can be used to tokenise the input stream line-wise:
 * \snippet test/snippet/range/view/take_line.cpp tokenise
 */
inline auto constexpr take_line = view::take_until_and_consume(is_char<'\r'> || is_char<'\n'>);

// ============================================================================
//  view::take_line_or_throw (adaptor instance definition)
// ============================================================================

/*!\brief A view adaptor that returns a single line from the underlying range (throws if there is no end-of-line).
 * \throws seqan3::unexpected_end_of_input If the underlying range contains no end-of-line marker.
 * \ingroup view
 *
 * \copydetails seqan3::view::take_line
 */
inline auto constexpr take_line_or_throw = view::take_until_or_throw_and_consume(is_char<'\r'> || is_char<'\n'>);

//!\}

} // namespace seqan3::view
