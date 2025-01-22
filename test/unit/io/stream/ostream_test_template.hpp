// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <concepts>
#include <fstream>
#include <iostream>
#include <string>

#include <seqan3/io/stream/concept.hpp>
#include <seqan3/test/tmp_directory.hpp>
#include <seqan3/test/zlib_skip.hpp>

template <typename T>
class ostream : public ::testing::Test
{};

inline std::string const uncompressed{"The quick brown fox jumps over the lazy dog"};

TYPED_TEST_SUITE_P(ostream);

TYPED_TEST_P(ostream, concept_check)
{
    EXPECT_TRUE((seqan3::output_stream_over<TypeParam, char>));
}

TYPED_TEST_P(ostream, output)
{
    seqan3::test::tmp_directory tmp{};
    auto filename = tmp.path() / "ostream_test";

    {
        std::ofstream of{filename};
        TypeParam ogzf{of};

        ogzf << uncompressed << std::flush;
    }

    std::ifstream fi{filename, std::ios::binary};
    std::string buffer{std::istreambuf_iterator<char>{fi}, std::istreambuf_iterator<char>{}};

    if constexpr (TestFixture::zero_out_os_byte)
    {
        buffer[9] = '\x00'; // zero-out the OS byte.
        SEQAN3_TEST_GTEST_SKIP_ZLIB_DEFLATE;
    }

    EXPECT_EQ(buffer, TestFixture::compressed);
}

TYPED_TEST_P(ostream, output_type_erased)
{
    seqan3::test::tmp_directory tmp{};
    auto filename = tmp.path() / "ostream_test";

    {
        std::ofstream of{filename};

        std::unique_ptr<std::ostream> ogzf{new TypeParam{of}};

        *ogzf << uncompressed << std::flush;
    }

    std::ifstream fi{filename, std::ios::binary};
    std::string buffer{std::istreambuf_iterator<char>{fi}, std::istreambuf_iterator<char>{}};

    if constexpr (TestFixture::zero_out_os_byte)
    {
        buffer[9] = '\x00'; // zero-out the OS byte.
        SEQAN3_TEST_GTEST_SKIP_ZLIB_DEFLATE;
    }

    EXPECT_EQ(buffer, TestFixture::compressed);
}

REGISTER_TYPED_TEST_SUITE_P(ostream, concept_check, output, output_type_erased);
