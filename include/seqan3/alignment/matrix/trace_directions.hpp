// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides the declaration of seqan3::detail::trace_directions.
 */

#pragma once

#include <seqan3/core/add_enum_bitwise_operators.hpp>
#include <seqan3/core/debug_stream.hpp>

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
    //!\brief Trace comes from the above entry.
    up        = 0b00010,
    //!\brief Trace comes from the left entry.
    left      = 0b00100,
    //!\brief Trace comes from the above entry, while opening the gap.
    up_open   = 0b01000,
    //!\brief Trace comes from the left entry, while opening the gap.
    left_open = 0b10000
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
 */
inline debug_stream_type & operator<<(debug_stream_type & s, detail::trace_directions const trace)
{
    static char const * unicode[8]{"↺", "↖", "↑", "↖↑", "←", "↖←", "↑←", "↖↑←"};
    static char const * csv[8]{"N", "D", "U", "DU", "L", "DL", "UL", "DUL"};

    bool is_unicode = (s.flags2() & fmtflags2::utf8) == fmtflags2::utf8;
    auto const & trace_dir = is_unicode ? unicode : csv;

    s << trace_dir[static_cast<size_t>(trace) % 8u];
    return s;
}

} // namespace seqan3
