// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <sstream>

#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/core/debug_stream.hpp>

using seqan3::operator|;

static constexpr seqan3::detail::trace_directions N = seqan3::detail::trace_directions::none;
static constexpr seqan3::detail::trace_directions D = seqan3::detail::trace_directions::diagonal;
static constexpr seqan3::detail::trace_directions u = seqan3::detail::trace_directions::up;
static constexpr seqan3::detail::trace_directions l = seqan3::detail::trace_directions::left;
static constexpr seqan3::detail::trace_directions U = seqan3::detail::trace_directions::up_open;
static constexpr seqan3::detail::trace_directions L = seqan3::detail::trace_directions::left_open;

TEST(debug_stream_test, ascii)
{
    std::stringstream s{};
    seqan3::debug_stream_type stream{s};
    stream << N << ";" << D << ";" << U << ";" << L << ";" << (D|U) << ";" << (D|L) << ";" << (U|L) << ";" << (D|U|L);
    stream << ";" << u << ";" << l << ";" << (D|u) << ";" << (D|u|l) << ";" << (D|U|u|L|l);

    EXPECT_EQ(s.str(), "N;D;U;L;DU;DL;UL;DUL;u;l;Du;Dul;DUuLl");
}

TEST(debug_stream_test, unicode)
{
    std::stringstream s{};
    seqan3::debug_stream_type stream{s};
    stream << seqan3::fmtflags2::utf8;
    stream << N << ";" << D << ";" << U << ";" << L << ";" << (D|U) << ";" << (D|L) << ";" << (U|L) << ";" << (D|U|L);
    stream << ";" << u << ";" << l << ";" << (D|u) << ";" << (D|u|l) << ";" << (D|U|u|L|l);

    EXPECT_EQ(s.str(), "↺;↖;↑;←;↖↑;↖←;↑←;↖↑←;⇡;⇠;↖⇡;↖⇡⇠;↖↑⇡←⇠");
}
