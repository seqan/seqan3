// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::rank_to.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/range/view/deep.hpp>
#include <seqan3/std/ranges>

namespace seqan3::view
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief                A view over an alphabet, given a range of ranks.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \tparam alphabet_t   The alphabet to convert to; must satisfy seqan3::WritableSemialphabet.
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of converted elements. See below for the properties of the returned range.
 * \ingroup view
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/view/rank_to.hpp>
 * ```
 *
 * ### View properties
 *
 * This view is a **deep view:** Given a range-of-range as input (as opposed to just a range), it will apply
 * the transformation on the innermost range (instead of the outermost range).
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |---------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                        |
 * | std::ranges::ForwardRange       |                                       | *preserved*                                        |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                                        |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                                        |
 * | std::ranges::ContiguousRange    |                                       | *lost*                                             |
 * |                                 |                                       |                                                    |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                       |
 * | std::ranges::View               |                                       | *guaranteed*                                       |
 * | std::ranges::SizedRange         |                                       | *preserved*                                        |
 * | std::ranges::CommonRange        |                                       | *preserved*                                        |
 * | std::ranges::OutputRange        |                                       | *lost*                                             |
 * | seqan3::ConstIterableRange      |                                       | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | seqan3::reference_t             | seqan3::alphabet_rank_t<alphabet_t>   | `alphabet_t`                                       |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ###Example
 * \snippet test/snippet/range/view/rank_char.cpp rank_to
 * \hideinitializer
 */
template <typename alphabet_type>
//!\cond
    requires WritableSemialphabet<alphabet_type>
//!\endcond
inline auto const rank_to = deep{std::view::transform(
[] (alphabet_rank_t<alphabet_type> const in) -> alphabet_type
{
    return assign_rank_to(in, alphabet_type{});
})};

//!\}

} // namespace seqan3::view
