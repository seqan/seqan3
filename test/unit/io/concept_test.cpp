// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <fstream>
#include <istream>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include <seqan3/io/stream/concept.hpp>

TEST(io, input_stream_over)
{
    EXPECT_TRUE((seqan3::input_stream_over<std::istream, char>));
    EXPECT_FALSE((seqan3::input_stream_over<std::ostream, char>));
    EXPECT_TRUE((seqan3::input_stream_over<std::iostream, char>));
    EXPECT_TRUE((seqan3::input_stream_over<std::ifstream, char>));
    EXPECT_FALSE((seqan3::input_stream_over<std::ofstream, char>));
    EXPECT_TRUE((seqan3::input_stream_over<std::fstream, char>));
    EXPECT_TRUE((seqan3::input_stream_over<std::istringstream, char>));
    EXPECT_FALSE((seqan3::input_stream_over<std::ostringstream, char>));
    EXPECT_TRUE((seqan3::input_stream_over<std::stringstream, char>));
    EXPECT_FALSE((seqan3::input_stream_over<std::vector<char>, char>));
    EXPECT_FALSE((seqan3::input_stream_over<std::vector<char>, char>));
    EXPECT_FALSE((seqan3::input_stream_over<std::vector<char>, char>));
    EXPECT_FALSE((seqan3::input_stream_over<std::string, char>));
    EXPECT_FALSE((seqan3::input_stream_over<std::string, char>));
    EXPECT_FALSE((seqan3::input_stream_over<std::string, char>));
}

TEST(io, output_stream_over)
{
    EXPECT_FALSE((seqan3::output_stream_over<std::istream, char>));
    EXPECT_TRUE((seqan3::output_stream_over<std::ostream, char>));
    EXPECT_TRUE((seqan3::output_stream_over<std::iostream, char>));
    EXPECT_FALSE((seqan3::output_stream_over<std::ifstream, char>));
    EXPECT_TRUE((seqan3::output_stream_over<std::ofstream, char>));
    EXPECT_TRUE((seqan3::output_stream_over<std::fstream, char>));
    EXPECT_FALSE((seqan3::output_stream_over<std::istringstream, char>));
    EXPECT_TRUE((seqan3::output_stream_over<std::ostringstream, char>));
    EXPECT_TRUE((seqan3::output_stream_over<std::stringstream, char>));
    EXPECT_FALSE((seqan3::output_stream_over<std::vector<char>, char>));
    EXPECT_FALSE((seqan3::output_stream_over<std::vector<char>, char>));
    EXPECT_FALSE((seqan3::output_stream_over<std::vector<char>, char>));
    EXPECT_FALSE((seqan3::output_stream_over<std::string, char>));
    EXPECT_FALSE((seqan3::output_stream_over<std::string, char>));
    EXPECT_FALSE((seqan3::output_stream_over<std::string, char>));
}
