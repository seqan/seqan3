// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/std/ranges>

namespace seqan3
{
/*!\name Formatted output overloads
 * \{
 */
} // namespace seqan3

namespace seqan3::detail
{

//!\brief Helper function to print elements of a tuple separately.
template <typename char_t, typename tuple_t, std::size_t ...I>
void print_tuple(debug_stream_type<char_t> & s, tuple_t && t, std::index_sequence<I...> const &)
{
    using std::get;
    s << '(';
    ((s << (I == 0 ? "" : ",") << get<I>(t)), ...);
    s << ')';
}

} // namespace seqan3::detail

namespace seqan3
{

/*!\brief All tuples can be printed by printing their elements separately.
 * \tparam tuple_t Type of the tuple to be printed; must model seqan3::tuple_like.
 * \param s The seqan3::debug_stream.
 * \param t The tuple.
 * \relates seqan3::debug_stream_type
 */
template <typename tuple_t, typename char_t>
//!\cond
    requires (!std::ranges::input_range<tuple_t>) &&
             (!alphabet<tuple_t>) && // exclude alphabet_tuple_base
             tuple_like<remove_cvref_t<tuple_t>>
//!\endcond
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, tuple_t && t)
{
    detail::print_tuple(s, std::forward<tuple_t>(t),
                        std::make_index_sequence<std::tuple_size_v<remove_cvref_t<tuple_t>>>{});
    return s;
}

//!\}

} // namespace seqan3
