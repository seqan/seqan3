// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/detail/default_simd_backend.hpp>
#include <seqan3/core/simd/detail/default_simd_length.hpp>

#include <iostream>
#include <type_traits>

TEST(default_simd_length, int8_t)
{
    using scalar_t = int8_t;
    constexpr size_t max_length = seqan3::detail::default_simd_max_length<seqan3::detail::builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 16 ); break;
        case 32: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 32 ); break;
        case 64: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 64 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}

TEST(default_simd_length, int16_t)
{
    using scalar_t = int16_t;
    constexpr size_t max_length = seqan3::detail::default_simd_max_length<seqan3::detail::builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 8 ); break;
        case 32: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 16 ); break;
        case 64: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 32 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}

TEST(default_simd_length, int32_t)
{
    using scalar_t = int32_t;
    constexpr size_t max_length = seqan3::detail::default_simd_max_length<seqan3::detail::builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 4 ); break;
        case 32: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 8 ); break;
        case 64: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 16 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}

TEST(default_simd_length, int64_t)
{
    using scalar_t = int64_t;
    constexpr size_t max_length = seqan3::detail::default_simd_max_length<seqan3::detail::builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 2 ); break;
        case 32: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 4 ); break;
        case 64: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 8 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}

TEST(default_simd_length, uint8_t)
{
    using scalar_t = uint8_t;
    constexpr size_t max_length = seqan3::detail::default_simd_max_length<seqan3::detail::builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 16 ); break;
        case 32: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 32 ); break;
        case 64: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 64 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}

TEST(default_simd_length, uint16_t)
{
    using scalar_t = uint16_t;
    constexpr size_t max_length = seqan3::detail::default_simd_max_length<seqan3::detail::builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 8 ); break;
        case 32: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 16 ); break;
        case 64: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 32 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}

TEST(default_simd_length, uint32_t)
{
    using scalar_t = uint32_t;
    constexpr size_t max_length = seqan3::detail::default_simd_max_length<seqan3::detail::builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 4 ); break;
        case 32: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 8 ); break;
        case 64: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 16 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}

TEST(default_simd_length, uint64_t)
{
    using scalar_t = uint64_t;
    constexpr size_t max_length = seqan3::detail::default_simd_max_length<seqan3::detail::builtin_simd>;

    switch (max_length)
    {
        case 0: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 1 ); break;
        case 16: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 2 ); break;
        case 32: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 4 ); break;
        case 64: EXPECT_EQ((seqan3::detail::default_simd_length<scalar_t, seqan3::detail::builtin_simd>), 8 ); break;
        default: FAIL() << "Unsupported max_length";
    }
}
