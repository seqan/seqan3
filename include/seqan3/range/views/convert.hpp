// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::views::convert.
 */

#pragma once

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3::views
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view that converts each element in the input range (implicitly or via `static_cast`).
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of converted elements. See below for the properties of the returned range.
 * \ingroup views
 *
 * \details
 *
 * \header_file{seqan3/range/views/convert.hpp}
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)  |
 * |----------------------------------|:-------------------------------------:|:-------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                     |
 * | std::ranges::forward_range       |                                       | *preserved*                     |
 * | std::ranges::bidirectional_range |                                       | *preserved*                     |
 * | std::ranges::random_access_range |                                       | *preserved*                     |
 * | std::ranges::contiguous_range    |                                       | *lost*¹                         |
 * |                                  |                                       |                                 |
 * | std::ranges::viewable_range      | *required*                            | *guaranteed*                    |
 * | std::ranges::view                |                                       | *guaranteed*                    |
 * | std::ranges::sized_range         |                                       | *preserved*                     |
 * | std::ranges::common_range        |                                       | *preserved*                     |
 * | std::ranges::output_range        |                                       | *lost*¹                         |
 * | seqan3::const_iterable_range     |                                       | *preserved*                     |
 * |                                  |                                       |                                 |
 * | std::ranges::range_reference_t   | seqan3::convertible_to<out_t>         | `out_t`                         |
 *
 * ¹ These are preserved if `out_t` is the same as `std::ranges::range_reference_t<urng_t>`, i.e. no conversion
 * takes place.
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * Convert from `int` to `bool`:
 * \include test/snippet/range/views/convert_int_to_bool.cpp
 *
 * Convert from seqan3::dna15 to seqan3::dna5:
 * \include test/snippet/range/views/convert_15_to_5.cpp
 * \hideinitializer
 */
template <typename out_t>
auto const convert = std::views::transform([] (auto && in) -> out_t
{
    static_assert(std::convertible_to<decltype(in) &&, out_t> || explicitly_convertible_to<decltype(in) &&, out_t>,
                  "The reference type of the input to views::convert is not convertible to out_t.");

    if constexpr (std::convertible_to<decltype(in) &&, out_t>)
        return std::forward<decltype(in)>(in);
    else
        return static_cast<out_t>(std::forward<decltype(in)>(in));
});

//!\}

} // namespace seqan3::views
