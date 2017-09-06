// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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
// ==========================================================================

#include <gtest/gtest.h>

#include <seqan3/simd/concept.hpp>
#include <seqan3/simd/detail/default_simd_backend.hpp>
#include <seqan3/simd/detail/default_simd_length.hpp>

#include <iostream>
#include <type_traits>

using namespace seqan3;
using namespace seqan3::detail;

TEST(default_simd_length, int8_t)
{
    using scalar_t = int8_t;
    constexpr size_t max_length = default_simd_max_length<builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 16 ); break;
        case 32: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 32 ); break;
        case 64: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 64 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}

TEST(default_simd_length, int16_t)
{
    using scalar_t = int16_t;
    constexpr size_t max_length = default_simd_max_length<builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 8 ); break;
        case 32: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 16 ); break;
        case 64: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 32 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}

TEST(default_simd_length, int32_t)
{
    using scalar_t = int32_t;
    constexpr size_t max_length = default_simd_max_length<builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 4 ); break;
        case 32: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 8 ); break;
        case 64: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 16 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}

TEST(default_simd_length, int64_t)
{
    using scalar_t = int64_t;
    constexpr size_t max_length = default_simd_max_length<builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 2 ); break;
        case 32: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 4 ); break;
        case 64: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 8 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}

TEST(default_simd_length, uint8_t)
{
    using scalar_t = uint8_t;
    constexpr size_t max_length = default_simd_max_length<builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 16 ); break;
        case 32: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 32 ); break;
        case 64: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 64 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}

TEST(default_simd_length, uint16_t)
{
    using scalar_t = uint16_t;
    constexpr size_t max_length = default_simd_max_length<builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 8 ); break;
        case 32: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 16 ); break;
        case 64: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 32 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}

TEST(default_simd_length, uint32_t)
{
    using scalar_t = uint32_t;
    constexpr size_t max_length = default_simd_max_length<builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 4 ); break;
        case 32: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 8 ); break;
        case 64: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 16 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}

TEST(default_simd_length, uint64_t)
{
    using scalar_t = uint64_t;
    constexpr size_t max_length = default_simd_max_length<builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 2 ); break;
        case 32: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 4 ); break;
        case 64: EXPECT_EQ((default_simd_length<scalar_t, builtin_simd>), 8 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}
