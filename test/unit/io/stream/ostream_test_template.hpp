// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/concepts>
#include <fstream>
#include <iostream>
#include <string>

#include <seqan3/io/stream/concept.hpp>
#include <seqan3/test/tmp_directory.hpp>

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
        buffer[9] = '\x00'; // zero-out the OS byte.

    EXPECT_EQ(buffer, TestFixture::compressed);
    tmp.clean();
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
        buffer[9] = '\x00'; // zero-out the OS byte.

    EXPECT_EQ(buffer, TestFixture::compressed);
    tmp.clean();
}

REGISTER_TYPED_TEST_SUITE_P(ostream, concept_check, output, output_type_erased);
