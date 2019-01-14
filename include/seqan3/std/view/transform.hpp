// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::transform.
 */

#pragma once

#include <range/v3/view/transform.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3::view
{

/*!\brief A range adaptor that takes a invocable and returns a view of the elements with the invocable applied.
 * \tparam           urng_t The type of the range being processed. See below for requirements. [template parameter
 *                          is omitted in pipe notation]
 * \tparam      invocable_t The type of the invocable, must satisfy std::Invocable.
 * \param[in]        urange The range being processed. [parameter is omitted in pipe notation]
 * \param[in,out] invocable The invocable (usually a lambda function).
 * \returns A range of the elements produced by applied the invocable to each element in the underlying range.
 * \ingroup core_view
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                       |
 * |---------------------------------|:-------------------------------------:|:----------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                          |
 * | std::ranges::ForwardRange       |                                       | *preserved*                                          |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                                          |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                                          |
 * | std::ranges::ContiguousRange    |                                       | *lost*                                               |
 * |                                 |                                       |                                                      |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                         |
 * | std::ranges::View               |                                       | *guaranteed*                                         |
 * | std::ranges::SizedRange         |                                       | *preserved*                                          |
 * | std::ranges::CommonRange        |                                       | *preserved*                                          |
 * | std::ranges::OutputRange        |                                       | *lost*                                               |
 * | seqan3::const_iterable_concept  |                                       | *preserved*                                          |
 * |                                 |                                       |                                                      |
 * | seqan3::reference_t             |                                       | `decltype(invocable(seqan3::reference_t<urng_t>{}))` |
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
 * auto v1 = s | view::transform([](dna4 const l) { return to_char(l); }); // == "ACTNTGATAN"
 * ```
 * \hideinitializer
 */

inline constexpr auto transform = ranges::view::transform;

} // namespace seqan3::view
