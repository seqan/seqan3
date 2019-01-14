// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Provides seqan3::view::get.
 */

#pragma once

#include <seqan3/std/view/transform.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/concept/tuple.hpp>

namespace seqan3::view
{
/*!\name General purpose views
 * \{
 */

/*!\brief               A view calling std::get on each element in a range.
 * \tparam size_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \tparam index        The index to get.
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of elements where every element is the result of calling std::get<index> on the underlying element.
                        See below for the properties of the returned range.
 * \ingroup view
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                          |
 * |---------------------------------|:-------------------------------------:|:-------------------------------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                                             |
 * | std::ranges::ForwardRange       |                                       | *preserved*                                             |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                                             |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                                             |
 * | std::ranges::ContiguousRange    |                                       | *lost*                                                  |
 * |                                 |                                       |                                                         |
 * | std::ranges::ViewableRange      | *required*                            | *preserved*                                             |
 * | std::ranges::View               |                                       | *preserved*                                             |
 * | std::ranges::SizedRange         |                                       | *preserved*                                             |
 * | std::ranges::CommonRange        |                                       | *preserved*                                             |
 * | std::ranges::OutputRange        |                                       | *preserved*                                             |
 * | seqan3::const_iterable_concept  |                                       | *preserved*                                             |
 * |                                 |                                       |                                                         |
 * | seqan3::reference_t             | seqan3::tuple_like_concept            | std::tuple_element_t<index, seqan3::reference_t<urng_t>>|
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * \snippet test/snippet/range/view/get.cpp usage
 * \hideinitializer
 */
template <size_t index>
inline auto const get = view::transform([] (auto && in) -> decltype(auto)
{
    using std::get;
    using seqan3::get;
    static_assert(tuple_like_concept<decltype(in)>,
                  "You may only pass ranges to view::get whose reference_t models the tuple_like_concept.");
    return get<index>(std::forward<decltype(in)>(in));
});

//!\}

} // namespace seqan3::view
