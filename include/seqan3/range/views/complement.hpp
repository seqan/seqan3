// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::views::complement.
 */

#pragma once

#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/range/views/deep.hpp>
#include <seqan3/std/ranges>

namespace seqan3::views
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief               A view that converts a range of nucleotides to their complement.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             A range of converted elements. See below for the properties of the returned range.
 * \ingroup views
 *
 * \details
 *
 * \header_file{seqan3/range/views/complement.hpp}
 *
 * Calls seqan3::nucleotide_alphabet::complement() on every element of the input range.
 *
 * ### View properties
 *
 * This view is a **deep view** Given a range-of-range as input (as opposed to just a range), it will apply
 * the transformation on the innermost range (instead of the outermost range).
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                                  |
 * |----------------------------------|:-------------------------------------:|:---------------------------------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                                                     |
 * | std::ranges::forward_range       |                                       | *preserved*                                                     |
 * | std::ranges::bidirectional_range |                                       | *preserved*                                                     |
 * | std::ranges::random_access_range |                                       | *preserved*                                                     |
 * | std::ranges::contiguous_range    |                                       | *lost*                                                          |
 * |                                  |                                       |                                                                 |
 * | std::ranges::viewable_range      | *required*                            | *guaranteed*                                                    |
 * | std::ranges::view                |                                       | *guaranteed*                                                    |
 * | std::ranges::sized_range         |                                       | *preserved*                                                     |
 * | std::ranges::common_range        |                                       | *preserved*                                                     |
 * | std::ranges::output_range        |                                       | *lost*                                                          |
 * | seqan3::const_iterable_range     |                                       | *preserved*                                                     |
 * |                                  |                                       |                                                                 |
 * | std::ranges::range_reference_t   | seqan3::nucleotide_alphabet           | std::remove_reference_t<std::ranges::range_reference_t<urng_t>> |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Example
 *
 * \include test/snippet/range/views/complement.cpp
 * \hideinitializer
 */

inline auto const complement = deep{std::views::transform([] (auto const in)
{
    static_assert(nucleotide_alphabet<decltype(in)>,
                  "The innermost value type must satisfy the nucleotide_alphabet.");
    // call element-wise complement from the nucleotide_alphabet
    return seqan3::complement(in);
})};

//!\}

} // namespace seqan3::views
