// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/debug_stream/debug_stream_type.hpp>
#include <seqan3/io/stream/concept.hpp>

namespace seqan3
{
/*!\name Formatted output overloads
 * \{
 */
/*!\brief All alphabets can be printed to the seqan3::debug_stream by their char representation.
 * \tparam alphabet_t Type of the alphabet to be printed; must model seqan3::alphabet.
 * \param s The seqan3::debug_stream.
 * \param l The alphabet letter.
 * \relates seqan3::debug_stream_type
 */
template <typename char_t, alphabet alphabet_t>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, alphabet_t && l)
    requires (!output_stream_over<std::basic_ostream<char_t>, alphabet_t>)
{
    return s << to_char(l);
}

// forward declare seqan3::mask
class mask;

/*!\brief Overload for the seqan3::mask alphabet.
 * \tparam char_t Type char type of the debug_stream.
 * \param s The seqan3::debug_stream.
 * \param l The mask alphabet letter.
 * \relates seqan3::debug_stream_type
 */
template <typename char_t, seqan3::semialphabet alphabet_t>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, alphabet_t && l)
    requires std::same_as<std::remove_cvref_t<alphabet_t>, mask>
{
    return s << (l == alphabet_t{} ? "UNMASKED" : "MASKED");
}

//!\}

} // namespace seqan3
