// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::convert.
 */

#pragma once

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3::view
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view that converts each element in the input range (implicitly or via `static_cast`).
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of converted elements. See below for the properties of the returned range.
 * \ingroup view
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/view/convert.hpp>
 * ```
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)  |
 * |---------------------------------|:-------------------------------------:|:-------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                     |
 * | std::ranges::ForwardRange       |                                       | *preserved*                     |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                     |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                     |
 * |                                 |                                       |                                 |
 * | std::ranges::View               |                                       | *guaranteed*                    |
 * | std::ranges::SizedRange         |                                       | *preserved*                     |
 * | std::ranges::CommonRange        |                                       | *preserved*                     |
 * | std::ranges::OutputRange        |                                       | *lost*                          |
 * | seqan3::ConstIterableRange      |                                       | *preserved*                     |
 * |                                 |                                       |                                 |
 * | seqan3::reference_t             | seqan3::ConvertibleTo<out_t>          | `out_t`                         |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * Convert from `int` to `bool`:
 * \snippet test/snippet/range/view/convert.cpp int_to_bool
 *
 * Convert from seqan3::dna15 to seqan3::dna5:
 * \snippet test/snippet/range/view/convert.cpp 15_to_5
 * \hideinitializer
 */
template <typename out_t>
auto const convert = std::view::transform([] (auto const & in) -> out_t
{
    if constexpr (ImplicitlyConvertibleTo<std::remove_reference_t<decltype(in)>, out_t>)
        return in;
    else
        return static_cast<out_t>(in);
});

//!\}

} // namespace seqan3::view
