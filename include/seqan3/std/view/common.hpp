// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::common.
 */

#pragma once

#include <range/v3/view/common.hpp>

#include <seqan3/std/iterator>

namespace seqan3::view
{

/*!\brief A range adaptor that makes any range satisfy std::ranges::CommonRange (at the expense of some performance).
 * \tparam           urng_t The type of the range being processed. See below for requirements. [template parameter
 *                          is omitted in pipe notation]
 * \param[in]        urange The range being processed. [parameter is omitted in pipe notation]
 * \returns A view of the underlying range that is common – even if the underlying range is not.
 * \ingroup core_view
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
 * | std::ranges::SizedRange         |                                       | *preserved*                                        |
 * | std::ranges::CommonRange        |                                       | *guarenteed*                                       |
 * | std::ranges::OutputRange        |                                       | *preserved*                                        |
 * | seqan3::const_iterable_concept  |                                       | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | seqan3::reference_t             |                                       | seqan3::reference_t<urng_t>                        |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### STD module
 *
 * This entity will likely be part of C++20 or C++23. It's API will track current proposals and not be stable
 * within SeqAn3 releases. It is implemented via the range-v3 library or the standard library (if available).
 * Should it become clear that it will not become part of a future standard, it will migrate to a regular SeqAn3
 * module.
 *
 * ### Example
 *
 * ```cpp
 * dna5_vector s{"ACTNTGATAN"_dna5};
 * auto v1 = s | view::filter([](dna5 const l) { return (l != 'N'_dna5); }); // == "ACTTGATA"_dna5
 *
 * // this won't work (as of C++17), because std::find expects begin and end to be of the same type:
 * // auto it = std::find(begin(v1), end(v1), 'G'_dna5);
 *
 * // this will:
 * auto v2 = v1 | view::common;
 * auto it = std::find(begin(v2), end(v2), 'G'_dna5);
 * ```
 * \hideinitializer
 */
inline constexpr auto common = ranges::view::common;

} // namespace seqan3::view
