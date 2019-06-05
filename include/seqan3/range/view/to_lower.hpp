// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Tobias Loka <LokaT AT rki.de>
 * \brief Provides seqan3::view::to_lower.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/char_operations/transform.hpp>
#include <seqan3/range/view/deep.hpp>
#include <seqan3/std/ranges>

namespace seqan3::view
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view that calls seqan3::to_lower() on each element in the input range.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of converted elements. See below for the properties of the returned range.
 * \ingroup view
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/view/to_lower.hpp>
 * ```
 *
 * ### View properties
 *
 * This view is a **deep view:** Given a range-of-range as input (as opposed to just a range), it will apply
 * the transformation on the innermost range (instead of the outermost range).
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type) | `rrng_t` (returned range type)                          |
 * |---------------------------------|:--------------------------------:|:-------------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                       | *preserved*                                             |
 * | std::ranges::ForwardRange       |                                  | *preserved*                                             |
 * | std::ranges::BidirectionalRange |                                  | *preserved*                                             |
 * | std::ranges::RandomAccessRange  |                                  | *preserved*                                             |
 * | std::ranges::ContiguousRange    |                                  | *lost*                                                  |
 * |                                 |                                  |                                                         |
 * | std::ranges::ViewableRange      | *required*                       | *guaranteed*                                            |
 * | std::ranges::View               |                                  | *guaranteed*                                            |
 * | std::ranges::SizedRange         |                                  | *preserved*                                             |
 * | std::ranges::CommonRange        |                                  | *preserved*                                             |
 * | std::ranges::OutputRange        |                                  | *lost*                                                  |
 * | seqan3::ConstIterableRange      |                                  | *preserved*                                             |
 * |                                 |                                  |                                                         |
 * | seqan3::reference_t             | seqan3::Char                     | seqan3::remove_reference_t<seqan3::reference_t<urngt_>> |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 * \snippet test/snippet/range/view/to_case.cpp to_lower
 * \hideinitializer
 */
inline auto const to_lower = deep{std::view::transform([] (auto const in) noexcept
{
    static_assert(Char<remove_cvref_t<decltype(in)>>,
                  "The value type of seqan3::view::to_lower must model seqan3::Char.");
    return seqan3::to_lower(in);
})};

//!\}

} // namespace seqan3::view
