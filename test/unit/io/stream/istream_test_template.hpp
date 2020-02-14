// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <string>

#include <seqan3/io/stream/concept.hpp>

#include <seqan3/test/tmp_filename.hpp>

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
    seqan3::test::tmp_filename filename{"istream_test"};

    {
        std::ofstream fi{filename.get_path()};

        fi << TestFixture::compressed;
    }

    std::ifstream fi{filename.get_path(), std::ios::binary};
    TypeParam comp{fi};
    std::string buffer{std::istreambuf_iterator<char>{comp}, std::istreambuf_iterator<char>{}};

    EXPECT_EQ(buffer, uncompressed);
}

TYPED_TEST_P(istream, input_type_erased)
{
    seqan3::test::tmp_filename filename{"istream_test"};

    {
        std::ofstream fi{filename.get_path()};

        fi << TestFixture::compressed;
    }

    std::ifstream fi{filename.get_path(), std::ios::binary};
    std::unique_ptr<std::istream> comp{new TypeParam{fi}};
    std::string buffer{std::istreambuf_iterator<char>{*comp}, std::istreambuf_iterator<char>{}};

    EXPECT_EQ(buffer, uncompressed);
}

REGISTER_TYPED_TEST_SUITE_P(istream, concept_check, input, input_type_erased);
