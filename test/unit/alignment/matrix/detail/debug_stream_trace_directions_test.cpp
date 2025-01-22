// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <sstream>

#include <seqan3/alignment/matrix/detail/trace_directions.hpp>
#include <seqan3/core/debug_stream/debug_stream_type.hpp>

using seqan3::operator|;

static constexpr seqan3::detail::trace_directions N = seqan3::detail::trace_directions::none;
static constexpr seqan3::detail::trace_directions D = seqan3::detail::trace_directions::diagonal;
static constexpr seqan3::detail::trace_directions u = seqan3::detail::trace_directions::up;
static constexpr seqan3::detail::trace_directions l = seqan3::detail::trace_directions::left;
static constexpr seqan3::detail::trace_directions U = seqan3::detail::trace_directions::carry_up_open;
static constexpr seqan3::detail::trace_directions L = seqan3::detail::trace_directions::carry_left_open;

TEST(debug_stream_test, ascii)
{
    std::stringstream s{};
    seqan3::debug_stream_type stream{s};
    stream << N << ";" << D << ";" << U << ";" << L << ";" << (D | U) << ";" << (D | L) << ";" << (U | L) << ";"
           << (D | U | L);
    stream << ";" << u << ";" << l << ";" << (D | u) << ";" << (D | u | l) << ";" << (D | U | u | L | l);

    EXPECT_EQ(s.str(), "N;D;U;L;DU;DL;UL;DUL;u;l;Du;Dul;DUuLl");
}

TEST(debug_stream_test, unicode)
{
    std::stringstream s{};
    seqan3::debug_stream_type stream{s};
    stream << seqan3::fmtflags2::utf8;
    stream << N << ";" << D << ";" << U << ";" << L << ";" << (D | U) << ";" << (D | L) << ";" << (U | L) << ";"
           << (D | U | L);
    stream << ";" << u << ";" << l << ";" << (D | u) << ";" << (D | u | l) << ";" << (D | U | u | L | l);

    EXPECT_EQ(s.str(), "↺;↖;↑;←;↖↑;↖←;↑←;↖↑←;⇡;⇠;↖⇡;↖⇡⇠;↖↑⇡←⇠");
}

TEST(trace_directions_printer_test, std_stream)
{
    std::stringstream s{};
    seqan3::trace_directions_printer<seqan3::detail::trace_directions> printer{};
    printer(s, N);
    s << ';';
    printer(s, D);
    s << ';';
    printer(s, L);
    s << ';';
    printer(s, U);
    s << ';';
    printer(s, D | U | u | L | l);

    EXPECT_EQ(s.str(), "N;D;L;U;DUuLl");
}
