// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::views::rank_to.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/utility/views/deep.hpp>

namespace seqan3::views
{
/*!\brief                A view over an alphabet, given a range of ranks.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \tparam alphabet_t   The alphabet to convert to; must satisfy seqan3::writable_semialphabet.
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of converted elements. See below for the properties of the returned range.
 * \ingroup alphabet_views
 *
 * \details
 *
 * \header_file{seqan3/alphabet/views/rank_to.hpp}
 *
 * ### View properties
 *
 * This view is a **deep view** Given a range-of-range as input (as opposed to just a range), it will apply
 * the transformation on the innermost range (instead of the outermost range).
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |----------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                                        |
 * | std::ranges::forward_range       |                                       | *preserved*                                        |
 * | std::ranges::bidirectional_range |                                       | *preserved*                                        |
 * | std::ranges::random_access_range |                                       | *preserved*                                        |
 * | std::ranges::contiguous_range    |                                       | *lost*                                             |
 * |                                  |                                       |                                                    |
 * | std::ranges::viewable_range      | *required*                            | *guaranteed*                                       |
 * | std::ranges::view                |                                       | *guaranteed*                                       |
 * | std::ranges::sized_range         |                                       | *preserved*                                        |
 * | std::ranges::common_range        |                                       | *preserved*                                        |
 * | std::ranges::output_range        |                                       | *lost*                                             |
 * | seqan3::const_iterable_range     |                                       | *preserved*                                        |
 * |                                  |                                       |                                                    |
 * | std::ranges::range_reference_t   | seqan3::alphabet_rank_t<alphabet_t>   | `alphabet_t`                                       |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ###Example
 * \include test/snippet/alphabet/views/rank_to.cpp
 * \hideinitializer
 *
 * \stableapi{Since version 3.1.}
 */
template <typename alphabet_type>
    requires writable_semialphabet<alphabet_type>
inline auto const rank_to = deep{std::views::transform(
    [](alphabet_rank_t<alphabet_type> const in) -> alphabet_type
    {
        return assign_rank_to(in, alphabet_type{});
    })};

} // namespace seqan3::views
