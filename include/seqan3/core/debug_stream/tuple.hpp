// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <concepts>
#include <ranges>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/debug_stream/debug_stream_type.hpp>
#include <seqan3/utility/tuple/concept.hpp>

namespace seqan3
{
/*!\name Formatted output overloads
 * \{
 */
} // namespace seqan3

namespace seqan3::detail
{

//!\brief Helper function to print elements of a tuple separately.
template <typename char_t, typename tuple_t, std::size_t... I>
void print_tuple(debug_stream_type<char_t> & s, tuple_t && t, std::index_sequence<I...> const &)
{
    using std::get;
    s << '(';
    ((s << (I == 0 ? "" : ",") << get<I>(t)), ...);
    s << ')';
}

/*!\interface seqan3::detail::debug_streamable_tuple <>
 * \brief A helper concept to avoid ambiguous overloads with the debug stream operator for alignments.
 * \ingroup core_debug_stream
 * \tparam tuple_t The tuple type to print to the seqan3::detail::debug_stream_type.
 *
 * \details
 *
 * This concept requires that the given type is a seqan3::tuple_like type but neither a std::ranges::input_range nor
 * an alphabet (see seqan3::alphabet_tuple_base).
 */
//!\cond
template <typename tuple_t>
concept debug_streamable_tuple =
    !std::ranges::input_range<tuple_t> && !alphabet<tuple_t> && // exclude alphabet_tuple_base
    tuple_like<std::remove_cvref_t<tuple_t>>;
//!\endcond
} // namespace seqan3::detail

namespace seqan3
{

/*!\brief All tuples can be printed by printing their elements separately.
 * \tparam tuple_t Type of the tuple to be printed; must model seqan3::tuple_like.
 * \param s The seqan3::debug_stream.
 * \param t The tuple.
 * \relates seqan3::debug_stream_type
 */
template <typename char_t, typename tuple_t>
    requires (detail::debug_streamable_tuple<tuple_t>)
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, tuple_t && t)
{
    detail::print_tuple(s,
                        std::forward<tuple_t>(t),
                        std::make_index_sequence<std::tuple_size_v<std::remove_cvref_t<tuple_t>>>{});
    return s;
}

//!\}

} // namespace seqan3
