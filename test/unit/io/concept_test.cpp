// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <istream>
#include <ostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <seqan3/io/stream/concept.hpp>

using namespace seqan3;

TEST(io, stream_REMOVEME)
{
    EXPECT_FALSE((stream_REMOVEME<std::istream, char>));
    EXPECT_FALSE((stream_REMOVEME<std::ostream, char>));
    EXPECT_TRUE((stream_REMOVEME<std::iostream, char>));
    EXPECT_FALSE((stream_REMOVEME<std::ifstream, char>));
    EXPECT_FALSE((stream_REMOVEME<std::ofstream, char>));
    EXPECT_TRUE((stream_REMOVEME<std::fstream, char>));
    EXPECT_FALSE((stream_REMOVEME<std::istringstream, char>));
    EXPECT_FALSE((stream_REMOVEME<std::ostringstream, char>));
    EXPECT_TRUE((stream_REMOVEME<std::stringstream, char>));
    EXPECT_FALSE((stream_REMOVEME<std::vector<char>, char>));
    EXPECT_FALSE((stream_REMOVEME<std::vector<char>, char>));
    EXPECT_FALSE((stream_REMOVEME<std::vector<char>, char>));
    EXPECT_FALSE((stream_REMOVEME<std::string, char>));
    EXPECT_FALSE((stream_REMOVEME<std::string, char>));
    EXPECT_FALSE((stream_REMOVEME<std::string, char>));
}

TEST(io, input_stream_over)
{
    EXPECT_TRUE((input_stream_over<std::istream, char>));
    EXPECT_FALSE((input_stream_over<std::ostream, char>));
    EXPECT_TRUE((input_stream_over<std::iostream, char>));
    EXPECT_TRUE((input_stream_over<std::ifstream, char>));
    EXPECT_FALSE((input_stream_over<std::ofstream, char>));
    EXPECT_TRUE((input_stream_over<std::fstream, char>));
    EXPECT_TRUE((input_stream_over<std::istringstream, char>));
    EXPECT_FALSE((input_stream_over<std::ostringstream, char>));
    EXPECT_TRUE((input_stream_over<std::stringstream, char>));
    EXPECT_FALSE((input_stream_over<std::vector<char>, char>));
    EXPECT_FALSE((input_stream_over<std::vector<char>, char>));
    EXPECT_FALSE((input_stream_over<std::vector<char>, char>));
    EXPECT_FALSE((input_stream_over<std::string, char>));
    EXPECT_FALSE((input_stream_over<std::string, char>));
    EXPECT_FALSE((input_stream_over<std::string, char>));
}

TEST(io, output_stream_over)
{
    EXPECT_FALSE((output_stream_over<std::istream, char>));
    EXPECT_TRUE((output_stream_over<std::ostream, char>));
    EXPECT_TRUE((output_stream_over<std::iostream, char>));
    EXPECT_FALSE((output_stream_over<std::ifstream, char>));
    EXPECT_TRUE((output_stream_over<std::ofstream, char>));
    EXPECT_TRUE((output_stream_over<std::fstream, char>));
    EXPECT_FALSE((output_stream_over<std::istringstream, char>));
    EXPECT_TRUE((output_stream_over<std::ostringstream, char>));
    EXPECT_TRUE((output_stream_over<std::stringstream, char>));
    EXPECT_FALSE((output_stream_over<std::vector<char>, char>));
    EXPECT_FALSE((output_stream_over<std::vector<char>, char>));
    EXPECT_FALSE((output_stream_over<std::vector<char>, char>));
    EXPECT_FALSE((output_stream_over<std::string, char>));
    EXPECT_FALSE((output_stream_over<std::string, char>));
    EXPECT_FALSE((output_stream_over<std::string, char>));
}
