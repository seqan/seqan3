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

TEST(io, OStream)
{
    EXPECT_FALSE((OStream<std::istream, char>));
    EXPECT_TRUE((OStream<std::ostream, char>));
    EXPECT_TRUE((OStream<std::iostream, char>));
    EXPECT_FALSE((OStream<std::ifstream, char>));
    EXPECT_TRUE((OStream<std::ofstream, char>));
    EXPECT_TRUE((OStream<std::fstream, char>));
    EXPECT_FALSE((OStream<std::istringstream, char>));
    EXPECT_TRUE((OStream<std::ostringstream, char>));
    EXPECT_TRUE((OStream<std::stringstream, char>));
    EXPECT_FALSE((OStream<std::vector<char>, char>));
    EXPECT_FALSE((OStream<std::vector<char>, char>));
    EXPECT_FALSE((OStream<std::vector<char>, char>));
    EXPECT_FALSE((OStream<std::string, char>));
    EXPECT_FALSE((OStream<std::string, char>));
    EXPECT_FALSE((OStream<std::string, char>));
}
