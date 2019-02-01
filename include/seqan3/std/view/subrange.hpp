// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::subrange.
 */

#pragma once

#include <range/v3/view/subrange.hpp>

#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

namespace seqan3::view
{

/*!\brief Create a view from a pair of iterator and sentinel.
 * \tparam   it_t Type of the iterator; must satisfy std::Iterator.
 * \tparam  sen_t Type of the sentinel; must satisfy std::Sentinel with it_t.
 * \param[in]  it The iterator on the underlying range.
 * \param[in] sen The sentinel on the underlying range
 * \returns  A view of the elements between it_t and sen_t.
 * \ingroup core_view
 *
 * ### View properties
 *
 * This view is **source-only**, it can only be at the beginning of a pipe of range transformations.
 *
 * | range concepts and reference_t  | `rrng_t` (returned range type)                     |
 * |---------------------------------|:--------------------------------------------------:|
 * | std::ranges::InputRange         | *preserved*                                        |
 * | std::ranges::ForwardRange       | *preserved*                                        |
 * | std::ranges::BidirectionalRange | *preserved*                                        |
 * | std::ranges::RandomAccessRange  | *preserved*                                        |
 * | std::ranges::ContiguousRange    | *preserved*                                        |
 * |                                 |                                                    |
 * | std::ranges::ViewableRange      | *guaranteed*                                       |
 * | std::ranges::View               | *guaranteed*                                       |
 * | std::ranges::SizedRange         | *preserved*                                        |
 * | std::ranges::CommonRange        | *preserved*                                        |
 * | std::ranges::OutputRange        | *preserved*                                        |
 * | seqan3::const_iterable_concept  | *preserved*                                        |
 * |                                 |                                                    |
 * | seqan3::reference_t             | seqan3::value_type_t<it_t>                         |
 *
 * Preservation in this table refers to the properties of the iterator/sentinel pair.
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
 * \snippet test/snippet/std/view/subrange.cpp example
 * \hideinitializer
 */
template <std::Iterator it_t, std::Sentinel<it_t> sen_t>
using subrange = ::ranges::subrange<it_t, sen_t>;
//TODO change to ranges::subrange once that has arrived

} // namespace seqan3::view
