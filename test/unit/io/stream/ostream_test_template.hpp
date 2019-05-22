// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <fstream>
#include <iostream>
#include <string>

#include <seqan3/contrib/stream/bgzf_ostream.hpp>
#include <seqan3/contrib/stream/gz_ostream.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/std/concepts>
#include <seqan3/test/tmp_filename.hpp>

using namespace seqan3;

template <typename T>
class ostream : public ::testing::Test
{};

inline std::string const uncompressed{"The quick brown fox jumps over the lazy dog"};

TYPED_TEST_CASE_P(ostream);

TYPED_TEST_P(ostream, concept_check)
{
    EXPECT_TRUE((OStream<TypeParam, char>));
}

TYPED_TEST_P(ostream, output)
{
    test::tmp_filename filename{"ostream_test"};

    {
        std::ofstream of{filename.get_path()};
        TypeParam ogzf{of};

        ogzf << uncompressed << std::flush;
    }

    std::ifstream fi{filename.get_path(), std::ios::binary};
    std::string buffer{std::istreambuf_iterator<char>{fi}, std::istreambuf_iterator<char>{}};

    if constexpr (std::Same<TypeParam, contrib::gz_ostream> || std::Same<TypeParam, contrib::bgzf_ostream>)
        buffer[9] = '\x00'; // zero-out the OS byte.

    EXPECT_EQ(buffer, TestFixture::compressed);
}

TYPED_TEST_P(ostream, output_type_erased)
{
    test::tmp_filename filename{"ostream_test"};

    {
        std::ofstream of{filename.get_path()};

        std::unique_ptr<std::ostream> ogzf{new TypeParam{of}};

        *ogzf << uncompressed << std::flush;
    }

    std::ifstream fi{filename.get_path(), std::ios::binary};
    std::string buffer{std::istreambuf_iterator<char>{fi}, std::istreambuf_iterator<char>{}};

    if constexpr (std::Same<TypeParam, contrib::gz_ostream> || std::Same<TypeParam, contrib::bgzf_ostream>)
        buffer[9] = '\x00'; // zero-out the OS byte.

    EXPECT_EQ(buffer, TestFixture::compressed);
}

REGISTER_TYPED_TEST_CASE_P(ostream, concept_check, output, output_type_erased);
