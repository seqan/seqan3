// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::all.
 */

#pragma once

#include <range/v3/view/all.hpp>

#include <seqan3/std/iterator>

namespace seqan3::view
{

/*!\brief A view that safely wraps a container (you will likely not need to use this unless defining a new view).
 * \tparam    urng_t The type of the range being processed. See below for requirements.
 * \param[in] urange The range being processed.
 * \returns A view over the elements of the underlying range.
 * \ingroup core_view
 *
 * ### View properties
 *
 * This view is **source-only**, it can only be at the beginning of a pipe of range transformations. The
 * "underlying range" refers to the mandatory parameter of this adaptor, not a previous range in a pipe.
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |---------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::InputRange         |                                       | *preserved*                                        |
 * | std::ranges::ForwardRange       |                                       | *preserved*                                        |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                                        |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                                        |
 * | std::ranges::ContiguousRange    |                                       | *preserved*                                        |
 * |                                 |                                       |                                                    |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                       |
 * | std::ranges::View               |                                       | *guaranteed*                                       |
 * | std::ranges::SizedRange         |                                       | *preserved*                                        |
 * | std::ranges::CommonRange        |                                       | *preserved*                                        |
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
 * dna4_vector s{"ACTTTGATAN"_dna4};
 * auto v = view::all(s); // the same as view::subrange{begin(s), end(s)}
 * ```
 * \hideinitializer
 */

inline constexpr auto all = ranges::view::all;

} // namespace seqan3::view

namespace seqan3::detail
{

//!\brief A transformation trait that returns for a given type `t` the return type of seqan3::view::all(t{}).
template <typename t>
using view_all_t = ranges::view::all_t<t>;

} // namespace seqan3::detail
