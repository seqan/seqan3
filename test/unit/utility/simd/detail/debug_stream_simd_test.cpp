// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>
#include <sstream>

#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/detail/debug_stream_simd.hpp>
#include <seqan3/utility/simd/simd.hpp>

TEST(debug_stream_test, simd_rvalue)
{
    using simd_type = seqan3::simd::simd_type_t<int16_t, 8>;

    std::stringstream strstream;
    seqan3::debug_stream_type stream{strstream};

    stream << seqan3::simd::iota<simd_type>(1);
    EXPECT_EQ(strstream.str(), "[1,2,3,4,5,6,7,8]");
}

TEST(debug_stream_test, simd_lvalue)
{
    using simd_type = seqan3::simd::simd_type_t<int16_t, 8>;

    std::stringstream strstream;
    seqan3::debug_stream_type stream{strstream};

    simd_type simd = seqan3::simd::iota<simd_type>(1);
    stream << simd;
    EXPECT_EQ(strstream.str(), "[1,2,3,4,5,6,7,8]");
}
