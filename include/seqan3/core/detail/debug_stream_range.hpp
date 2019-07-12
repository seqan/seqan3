// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{
/*!\name Formatted output overloads
 * \{
 */
/*!\brief All input ranges can be printed to the seqan3::debug_stream element-wise (if their elements are printable).
 * \tparam rng_t Type of the range to be printed; must model std::ranges::InputRange.
 * \param s The seqan3::debug_stream.
 * \param r The input range.
 * \relates seqan3::debug_stream_type
 *
 * \details
 *
 * If the element type models seqan3::Alphabet (and is not an unsigned integer), the range is printed
 * just as if it were a string, i.e. <tt>std::vector<dna4>{'C'_dna4, 'G'_dna4, 'A'_dna4}</tt> is printed as "CGA".
 *
 * In all other cases the elements are comma separated and the range is enclosed in brackets, i.e.
 * `std::vector<int>{3, 1, 33, 7}` is printed as "[3,1,33,7]".
 */
template <std::ranges::InputRange rng_t, typename char_t>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, rng_t && r)
//!\cond
    requires !std::Same<remove_cvref_t<reference_t<rng_t>>, remove_cvref_t<rng_t>> && // prevent recursive instantiation
             requires (reference_t<rng_t> l) { { s << l }; } &&
             // exclude null-terminated strings:
             !(std::is_pointer_v<std::decay_t<rng_t>> &&
               std::Same<remove_cvref_t<reference_t<rng_t>>, char>)
//!\endcond
{
    if constexpr (Alphabet<reference_t<rng_t>> &&
                  !detail::is_uint_adaptation_v<remove_cvref_t<reference_t<rng_t>>>)
    {
        for (auto && l : r)
            s << l;
    }
    else
    {
        s << '[';
        auto b = begin(r);
        auto e = end(r);
        if (b != e)
        {
            s << *b;
            ++b;
        }
        while (b != e)
        {
            s << ',';
            s << *b;
            ++b;
        }
        s << ']';
    }

    return s;
}

//!\}

} // namespace seqan3
