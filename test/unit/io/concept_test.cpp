// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
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

TEST(io, IStream)
{
    EXPECT_TRUE((IStream<std::istream, char>));
    EXPECT_FALSE((IStream<std::ostream, char>));
    EXPECT_TRUE((IStream<std::iostream, char>));
    EXPECT_TRUE((IStream<std::ifstream, char>));
    EXPECT_FALSE((IStream<std::ofstream, char>));
    EXPECT_TRUE((IStream<std::fstream, char>));
    EXPECT_TRUE((IStream<std::istringstream, char>));
    EXPECT_FALSE((IStream<std::ostringstream, char>));
    EXPECT_TRUE((IStream<std::stringstream, char>));
    EXPECT_FALSE((IStream<std::vector<char>, char>));
    EXPECT_FALSE((IStream<std::vector<char>, char>));
    EXPECT_FALSE((IStream<std::vector<char>, char>));
    EXPECT_FALSE((IStream<std::string, char>));
    EXPECT_FALSE((IStream<std::string, char>));
    EXPECT_FALSE((IStream<std::string, char>));
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
