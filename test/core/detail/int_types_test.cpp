// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

#include <gtest/gtest.h>
#include <seqan3/core/detail/int_types.hpp>

using namespace seqan3;

TEST(int_types_test, min_viable_uint_t)
{
    using bool_1_t = detail::min_viable_uint_t<0ull>;
    using bool_2_t = detail::min_viable_uint_t<1ull>;
    using uint8_1_t = detail::min_viable_uint_t<2ull>;
    using uint8_2_t = detail::min_viable_uint_t<0xFFull>;
    using uint16_1_t = detail::min_viable_uint_t<0x100ull>;
    using uint16_2_t = detail::min_viable_uint_t<0xFFFFull>;
    using uint32_1_t = detail::min_viable_uint_t<0x10000ull>;
    using uint32_2_t = detail::min_viable_uint_t<0xFFFFFFFFull>;
    using uint64_1_t = detail::min_viable_uint_t<0x100000000ull>;
    using uint64_2_t = detail::min_viable_uint_t<0xFFFFFFFFFFFFFFFFull>;

    EXPECT_TRUE((std::is_same_v<bool_1_t, bool>));
    EXPECT_TRUE((std::is_same_v<bool_2_t, bool>));
    EXPECT_TRUE((std::is_same_v<uint8_1_t, uint8_t>));
    EXPECT_TRUE((std::is_same_v<uint8_2_t, uint8_t>));
    EXPECT_TRUE((std::is_same_v<uint16_1_t, uint16_t>));
    EXPECT_TRUE((std::is_same_v<uint16_2_t, uint16_t>));
    EXPECT_TRUE((std::is_same_v<uint32_1_t, uint32_t>));
    EXPECT_TRUE((std::is_same_v<uint32_2_t, uint32_t>));
    EXPECT_TRUE((std::is_same_v<uint64_1_t, uint64_t>));
    EXPECT_TRUE((std::is_same_v<uint64_2_t, uint64_t>));
}

TEST(int_types_test, min_viable_uint_v)
{
    auto bool_1_v = detail::min_viable_uint_v<0ull>;
    auto bool_2_v = detail::min_viable_uint_v<1ull>;
    auto uint8_1_v = detail::min_viable_uint_v<2ull>;
    auto uint8_2_v = detail::min_viable_uint_v<0xFFull>;
    auto uint16_1_v = detail::min_viable_uint_v<0x100ull>;
    auto uint16_2_v = detail::min_viable_uint_v<0xFFFFull>;
    auto uint32_1_v = detail::min_viable_uint_v<0x10000ull>;
    auto uint32_2_v = detail::min_viable_uint_v<0xFFFFFFFFull>;
    auto uint64_1_v = detail::min_viable_uint_v<0x100000000ull>;
    auto uint64_2_v = detail::min_viable_uint_v<0xFFFFFFFFFFFFFFFFull>;

    EXPECT_EQ(static_cast<uint64_t>(bool_1_v), 0ull);
    EXPECT_EQ(static_cast<uint64_t>(bool_2_v), 1ull);
    EXPECT_EQ(static_cast<uint64_t>(uint8_1_v), 2ull);
    EXPECT_EQ(static_cast<uint64_t>(uint8_2_v), 0xFFull);
    EXPECT_EQ(static_cast<uint64_t>(uint16_1_v), 0x100ull);
    EXPECT_EQ(static_cast<uint64_t>(uint16_2_v), 0xFFFFull);
    EXPECT_EQ(static_cast<uint64_t>(uint32_1_v), 0x10000ull);
    EXPECT_EQ(static_cast<uint64_t>(uint32_2_v), 0xFFFFFFFFull);
    EXPECT_EQ(static_cast<uint64_t>(uint64_1_v), 0x100000000ull);
    EXPECT_EQ(static_cast<uint64_t>(uint64_2_v), 0xFFFFFFFFFFFFFFFFull);

    EXPECT_TRUE((std::is_same_v<decltype(bool_1_v), bool>));
    EXPECT_TRUE((std::is_same_v<decltype(bool_2_v), bool>));
    EXPECT_TRUE((std::is_same_v<decltype(uint8_1_v), uint8_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint8_2_v), uint8_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint16_1_v), uint16_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint16_2_v), uint16_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint32_1_v), uint32_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint32_2_v), uint32_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint64_1_v), uint64_t>));
    EXPECT_TRUE((std::is_same_v<decltype(uint64_2_v), uint64_t>));
}
