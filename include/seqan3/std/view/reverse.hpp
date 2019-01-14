// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::reverse.
 */

#pragma once

#include <range/v3/view/reverse.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3::view
{

/*!\brief A range adaptor that presents the underlying range in reverse order.
 * \tparam           urng_t The type of the range being processed. See below for requirements. [template parameter
 *                          is omitted in pipe notation]
 * \param[in]        urange The range being processed. [parameter is omitted in pipe notation]
 * \returns A view of the elements of the underlying range in reverse order.
 * \ingroup core_view
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                       |
 * |---------------------------------|:-------------------------------------:|:----------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                          |
 * | std::ranges::ForwardRange       | *required*                            | *preserved*                                          |
 * | std::ranges::BidirectionalRange | *required*                            | *preserved*                                          |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                                          |
 * | std::ranges::ContiguousRange    |                                       | *lost*                                               |
 * |                                 |                                       |                                                      |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                         |
 * | std::ranges::View               |                                       | *guaranteed*                                         |
 * | std::ranges::SizedRange         |                                       | *preserved*                                          |
 * | std::ranges::CommonRange        |                                       | *preserved*                                          |
 * | std::ranges::OutputRange        |                                       | *lost*                                               |
 * | seqan3::const_iterable_concept  |                                       | *preserved if std::ranges::CommonRange<urng_t>*  |
 * |                                 |                                       |                                                      |
 * | seqan3::reference_t             |                                       | seqan3::reference_t<urng_t>                          |
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
 * std::string s{"ACTNTGATAN"};
 * auto v1 = s | view::reverse; // == "NATAGTNTCA"
 * ```
 * \hideinitializer
 */

inline constexpr auto reverse = ranges::view::reverse;

} // namespace seqan3::view
