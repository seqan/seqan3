// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides the declaration of seqan3::detail::trace_directions.
 */

#pragma once

#include <array>
#include <string_view>

#include <seqan3/core/add_enum_bitwise_operators.hpp>
#include <seqan3/core/debug_stream/debug_stream_type.hpp>
#include <seqan3/core/detail/template_inspection.hpp>

namespace seqan3::detail
{

/*!\brief The possible directions a trace can have. The values can be combined by the logical `|`-operator.
 * \ingroup alignment_matrix
 * \implements seqan3::enum_bitwise_operators
 * \implements seqan3::enum_bitwise_operators
 * \sa seqan3::enum_bitwise_operators enables combining enum values.
 * \sa seqan3::detail::alignment_trace_matrix implementations use this enum as matrix entry type.
 */
enum struct trace_directions : uint8_t
{
    //!\brief No trace
    none = 0b0'0000,
    //!\brief Trace comes from the diagonal entry.
    diagonal = 0b0'0001,
    //!\brief Trace comes from the above entry, while opening the gap.
    up_open = 0b0'0110,
    //!\brief Trace comes from the above entry.
    up = 0b0'0100,
    //!\brief Trace comes from the left entry, while opening the gap.
    left_open = 0b1'1000,
    //!\brief Trace comes from the left entry.
    left = 0b1'0000,
    //!\brief Carry bit for the last up open even if it is not the maximum value.
    carry_up_open = 0b0'0010,
    //!\brief Carry bit for the last left open even if it is not the maximum value.
    carry_left_open = 0b0'1000
};

} // namespace seqan3::detail

namespace seqan3
{
//!\cond DEV
//!\brief Enable bitwise operators for enum seqan3::detail::trace_directions.
//!\ingroup alignment_matrix
//!\sa seqan3::enum_bitwise_operators enables combining enum values.
template <>
inline constexpr bool add_enum_bitwise_operators<seqan3::detail::trace_directions> = true;
//!\endcond

/*!\brief Prints `trace_directions` as ascii or as utf8 to output stream.
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
 *
 * \ingroup alignment_matrix
 */
template <>
struct trace_directions_printer<detail::trace_directions>
{
private:
    //!\brief The unicode representation of the trace directions.
    static constexpr std::array<std::string_view, 32> unicode{
        "↺", "↖",  "↑",  "↖↑",  "⇡",  "↖⇡",  "↑⇡",  "↖↑⇡",  "←",  "↖←",  "↑←",  "↖↑←",  "⇡←",  "↖⇡←",  "↑⇡←",  "↖↑⇡←",
        "⇠", "↖⇠", "↑⇠", "↖↑⇠", "⇡⇠", "↖⇡⇠", "↑⇡⇠", "↖↑⇡⇠", "←⇠", "↖←⇠", "↑←⇠", "↖↑←⇠", "⇡←⇠", "↖⇡←⇠", "↑⇡←⇠", "↖↑⇡←⇠"};

    //!\brief The ascii representation of the trace directions.
    static constexpr std::array<std::string_view, 32> csv{
        "N", "D",  "U",  "DU",  "u",  "Du",  "Uu",  "DUu",  "L",  "DL",  "UL",  "DUL",  "uL",  "DuL",  "UuL",  "DUuL",
        "l", "Dl", "Ul", "DUl", "ul", "Dul", "Uul", "DUul", "Ll", "DLl", "ULl", "DULl", "uLl", "DuLl", "UuLl", "DUuLl"};

public:
    /*!\brief Prints the trace directions into the given stream.
     *
     * This overload is only available if the stream has a member function `flags2` that returns a `fmtflags2`.
     * Using the flags2() member function allows to print the trace with unicode characters if seqan3::fmtflags2::utf8
     * is set to the seqan3::debug_stream.
     *
     * \tparam stream_t The type of the stream.
     * \param stream The stream to print to.
     * \param trace The trace directions to print.
     */
    template <typename stream_t>
        requires detail::is_type_specialisation_of_v<stream_t, debug_stream_type>
    constexpr void operator()(stream_t & stream, detail::trace_directions const trace) const
    {
        print_impl(stream, stream.flags2(), trace);
    }

    /*!\brief Prints the trace directions into the given stream.
     *
     * This overload is only available if the stream has no member function `flags2`. In this case it will use
     * ascii characters to print the trace.
     *
     * \tparam stream_t The type of the stream.
     * \param stream The stream to print to.
     * \param trace The trace directions to print.
     */
    template <typename stream_t>
    constexpr void operator()(stream_t & stream, detail::trace_directions const trace) const
    {
        print_impl(stream, fmtflags2::none, trace);
    }

private:
    /*!\brief Prints the trace directions
     * \tparam stream_t The type of the stream.
     * \param stream The stream to print to.
     * \param flag The flags of the stream.
     * \param trace The trace directions to print.
     */
    template <typename stream_t>
    constexpr void print_impl(stream_t & stream, fmtflags2 const flag, detail::trace_directions const trace) const
    {
        bool const is_unicode = (flag & fmtflags2::utf8) == fmtflags2::utf8;
        auto const & trace_dir = is_unicode ? unicode : csv;

        stream << trace_dir[static_cast<size_t>(trace)];
    }
};

} // namespace seqan3
