// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <string>

#include <seqan3/io/stream/concept.hpp>
#include <seqan3/test/tmp_directory.hpp>

template <typename T>
class istream : public ::testing::Test
{};

inline std::string const uncompressed{"The quick brown fox jumps over the lazy dog"};

TYPED_TEST_SUITE_P(istream);

TYPED_TEST_P(istream, concept_check)
{
    EXPECT_TRUE((seqan3::input_stream_over<TypeParam, char>));
}

TYPED_TEST_P(istream, input)
{
    seqan3::test::tmp_directory tmp{};
    auto filename = tmp.path() / "istream_test";

    {
        std::ofstream fi{filename};

        fi << TestFixture::compressed;
    }

    std::ifstream fi{filename, std::ios::binary};
    TypeParam comp{fi};
    std::string buffer{std::istreambuf_iterator<char>{comp}, std::istreambuf_iterator<char>{}};

    EXPECT_EQ(buffer, uncompressed);
}

TYPED_TEST_P(istream, input_type_erased)
{
    seqan3::test::tmp_directory tmp{};
    auto filename = tmp.path() / "istream_test";

    {
        std::ofstream fi{filename};

        fi << TestFixture::compressed;
    }

    std::ifstream fi{filename, std::ios::binary};
    std::unique_ptr<std::istream> comp{new TypeParam{fi}};
    std::string buffer{std::istreambuf_iterator<char>{*comp}, std::istreambuf_iterator<char>{}};

    EXPECT_EQ(buffer, uncompressed);
}

REGISTER_TYPED_TEST_SUITE_P(istream, concept_check, input, input_type_erased);
