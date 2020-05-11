// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------


#include <gtest/gtest.h>

#include <seqan3/core/simd/debug_stream_simd.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>

#include <iostream>
#include <sstream>

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
