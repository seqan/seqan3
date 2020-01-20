// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides the declaration of seqan3::detail::trace_directions.
 */

#pragma once

#include <seqan3/core/add_enum_bitwise_operators.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>

namespace seqan3::detail
{

/*!\brief The possible directions a trace can have. The values can be combined by the logical `|`-operator.
 * \ingroup alignment_matrix
 * \sa seqan3::add_enum_bitwise_operators <seqan3::detail::trace_directions> enables combining enum values.
 * \sa seqan3::detail::alignment_trace_matrix implementations use this enum as matrix entry type.
 */
enum struct trace_directions : uint8_t
{
    //!\brief No trace
    none      = 0b00000,
    //!\brief Trace comes from the diagonal entry.
    diagonal  = 0b00001,
    //!\brief Trace comes from the above entry, while opening the gap.
    up_open   = 0b00010,
    //!\brief Trace comes from the above entry.
    up        = 0b00100,
    //!\brief Trace comes from the left entry, while opening the gap.
    left_open = 0b01000,
    //!\brief Trace comes from the left entry.
    left      = 0b10000
};

} // namespace seqan3::detail

namespace seqan3
{
//!\brief Enable bitwise operators for enum seqan3::detail::trace_directions.
//!\ingroup alignment_matrix
template <>
constexpr bool add_enum_bitwise_operators<seqan3::detail::trace_directions> = true;
} // namespace seqan3

namespace seqan3
{

/*!\brief All trace_directions can be printed as ascii or as utf8 to the seqan3::debug_stream.
 * \param s The seqan3::debug_stream.
 * \param trace The trace direction.
 * \relates seqan3::debug_stream_type
 *
 * \details
 *
 * The following table shows the printed symbol of a particular seqan3::detail::trace_directions:
 *
 * | trace direction                             | utf8 | ascii |
 * |:-------------------------------------------:|:----:|:-----:|
 * | seqan3::detail::trace_directions::none      | ↺    | N     |
 * | seqan3::detail::trace_directions::diagonal  | ↖    | D     |
 * | seqan3::detail::trace_directions::up_open   | ↑    | U     |
 * | seqan3::detail::trace_directions::up        | ⇡    | u     |
 * | seqan3::detail::trace_directions::left_open | ←    | L     |
 * | seqan3::detail::trace_directions::left      | ⇠    | l     |
 */
template <typename char_t>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, detail::trace_directions const trace)
{
    static char const * unicode[32]{ "↺",   "↖",   "↑",   "↖↑",   "⇡",   "↖⇡",   "↑⇡",   "↖↑⇡",
                                     "←",  "↖←",  "↑←",  "↖↑←",  "⇡←",  "↖⇡←",  "↑⇡←",  "↖↑⇡←",
                                     "⇠",  "↖⇠",  "↑⇠",  "↖↑⇠",  "⇡⇠",  "↖⇡⇠",  "↑⇡⇠",  "↖↑⇡⇠",
                                    "←⇠", "↖←⇠", "↑←⇠", "↖↑←⇠", "⇡←⇠", "↖⇡←⇠", "↑⇡←⇠", "↖↑⇡←⇠"};

    static char const * csv[32]{ "N",   "D",   "U",   "DU",   "u",   "Du",   "Uu",   "DUu",
                                 "L",  "DL",  "UL",  "DUL",  "uL",  "DuL",  "UuL",  "DUuL",
                                 "l",  "Dl",  "Ul",  "DUl",  "ul",  "Dul",  "Uul",  "DUul",
                                "Ll", "DLl", "ULl", "DULl", "uLl", "DuLl", "UuLl", "DUuLl"};

    bool is_unicode = (s.flags2() & fmtflags2::utf8) == fmtflags2::utf8;
    auto const & trace_dir = is_unicode ? unicode : csv;

    s << trace_dir[static_cast<size_t>(trace)];
    return s;
}

} // namespace seqan3
