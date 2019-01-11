// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
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

TEST(io, Stream)
{
    EXPECT_FALSE((Stream<std::istream, char>));
    EXPECT_FALSE((Stream<std::ostream, char>));
    EXPECT_TRUE((Stream<std::iostream, char>));
    EXPECT_FALSE((Stream<std::ifstream, char>));
    EXPECT_FALSE((Stream<std::ofstream, char>));
    EXPECT_TRUE((Stream<std::fstream, char>));
    EXPECT_FALSE((Stream<std::istringstream, char>));
    EXPECT_FALSE((Stream<std::ostringstream, char>));
    EXPECT_TRUE((Stream<std::stringstream, char>));
    EXPECT_FALSE((Stream<std::vector<char>, char>));
    EXPECT_FALSE((Stream<std::vector<char>, char>));
    EXPECT_FALSE((Stream<std::vector<char>, char>));
    EXPECT_FALSE((Stream<std::string, char>));
    EXPECT_FALSE((Stream<std::string, char>));
    EXPECT_FALSE((Stream<std::string, char>));
}

TEST(io, istream_concept)
{
    EXPECT_TRUE((istream_concept<std::istream, char>));
    EXPECT_FALSE((istream_concept<std::ostream, char>));
    EXPECT_TRUE((istream_concept<std::iostream, char>));
    EXPECT_TRUE((istream_concept<std::ifstream, char>));
    EXPECT_FALSE((istream_concept<std::ofstream, char>));
    EXPECT_TRUE((istream_concept<std::fstream, char>));
    EXPECT_TRUE((istream_concept<std::istringstream, char>));
    EXPECT_FALSE((istream_concept<std::ostringstream, char>));
    EXPECT_TRUE((istream_concept<std::stringstream, char>));
    EXPECT_FALSE((istream_concept<std::vector<char>, char>));
    EXPECT_FALSE((istream_concept<std::vector<char>, char>));
    EXPECT_FALSE((istream_concept<std::vector<char>, char>));
    EXPECT_FALSE((istream_concept<std::string, char>));
    EXPECT_FALSE((istream_concept<std::string, char>));
    EXPECT_FALSE((istream_concept<std::string, char>));
}

TEST(io, ostream_concept)
{
    EXPECT_FALSE((ostream_concept<std::istream, char>));
    EXPECT_TRUE((ostream_concept<std::ostream, char>));
    EXPECT_TRUE((ostream_concept<std::iostream, char>));
    EXPECT_FALSE((ostream_concept<std::ifstream, char>));
    EXPECT_TRUE((ostream_concept<std::ofstream, char>));
    EXPECT_TRUE((ostream_concept<std::fstream, char>));
    EXPECT_FALSE((ostream_concept<std::istringstream, char>));
    EXPECT_TRUE((ostream_concept<std::ostringstream, char>));
    EXPECT_TRUE((ostream_concept<std::stringstream, char>));
    EXPECT_FALSE((ostream_concept<std::vector<char>, char>));
    EXPECT_FALSE((ostream_concept<std::vector<char>, char>));
    EXPECT_FALSE((ostream_concept<std::vector<char>, char>));
    EXPECT_FALSE((ostream_concept<std::string, char>));
    EXPECT_FALSE((ostream_concept<std::string, char>));
    EXPECT_FALSE((ostream_concept<std::string, char>));
}
