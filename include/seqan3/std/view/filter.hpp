// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::filter.
 */

#pragma once

#include <range/v3/view/filter.hpp>

#include <seqan3/std/iterator>

namespace seqan3::view
{

/*!\brief A range adaptor that takes a predicate and returns a view of the elements that satisfy the predicate.
 * \tparam           urng_t The type of the range being processed. See below for requirements. [template parameter
 *                          is omitted in pipe notation]
 * \tparam      predicate_t The type of the predicate, must satisfy seqan3::predicate_concept.
 * \param[in]        urange The range being processed. [parameter is omitted in pipe notation]
 * \param[in,out] predicate The predicate.
 * \returns A range of those elements in the underlying range that satisfy the predicate.
 * \ingroup core_view
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
 * |---------------------------------|:-------------------------------------:|:--------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                        |
 * | std::ranges::ForwardRange       |                                       | *preserved*                                        |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                                        |
 * | std::ranges::RandomAccessRange  |                                       | *lost*                                             |
 * | std::ranges::ContiguousRange    |                                       | *lost*                                             |
 * |                                 |                                       |                                                    |
 * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                       |
 * | std::ranges::View               |                                       | *guaranteed*                                       |
 * | std::ranges::SizedRange         |                                       | *lost*                                             |
 * | std::ranges::CommonRange        |                                       | *lost*                                             |
 * | std::ranges::OutputRange        |                                       | *preserved*                                        |
 * | seqan3::const_iterable_concept  |                                       | *lost*                                             |
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
 * ```
 * \hideinitializer
 */

inline constexpr auto filter = ranges::view::filter;

} // namespace seqan3::view
