// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#include <gtest/gtest.h>

#include <vector>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <string>

#include <seqan3/range/container/all.hpp>

using namespace seqan3;

TEST(container_concept, container_concept)
{
    EXPECT_TRUE((seqan3::container_concept<std::array<char, 2>>));
    EXPECT_TRUE((seqan3::container_concept<std::list<char>>));
    EXPECT_FALSE((seqan3::container_concept<std::forward_list<char>>)); // `.size()` missing
    EXPECT_TRUE((seqan3::container_concept<std::vector<char>>));
    EXPECT_TRUE((seqan3::container_concept<std::deque<char>>));
    EXPECT_TRUE((seqan3::container_concept<std::string>));
}

TEST(container_concept, sequence_concept)
{
    EXPECT_FALSE((seqan3::sequence_concept<std::array<char, 2>>));
    EXPECT_TRUE((seqan3::sequence_concept<std::list<char>>));
    EXPECT_FALSE((seqan3::sequence_concept<std::forward_list<char>>));
    EXPECT_TRUE((seqan3::sequence_concept<std::vector<char>>));
    EXPECT_TRUE((seqan3::sequence_concept<std::deque<char>>));
    EXPECT_TRUE((seqan3::sequence_concept<std::string>));
}

TEST(container_concept, random_access_sequence_concept)
{
    EXPECT_FALSE((seqan3::random_access_sequence_concept<std::array<char, 2>>));
    EXPECT_FALSE((seqan3::random_access_sequence_concept<std::list<char>>));
    EXPECT_FALSE((seqan3::random_access_sequence_concept<std::forward_list<char>>));
    EXPECT_TRUE((seqan3::random_access_sequence_concept<std::vector<char>>));
    EXPECT_TRUE((seqan3::random_access_sequence_concept<std::deque<char>>));
    EXPECT_TRUE((seqan3::random_access_sequence_concept<std::string>));
}

TEST(container_concept, reservable_sequence_concept)
{
    EXPECT_FALSE((seqan3::reservable_sequence_concept<std::array<char, 2>>));
    EXPECT_FALSE((seqan3::reservable_sequence_concept<std::list<char>>));
    EXPECT_FALSE((seqan3::reservable_sequence_concept<std::forward_list<char>>));
    EXPECT_TRUE((seqan3::reservable_sequence_concept<std::vector<char>>));
    EXPECT_FALSE((seqan3::reservable_sequence_concept<std::deque<char>>));
    EXPECT_TRUE((seqan3::reservable_sequence_concept<std::string>));
}

/* Check the SDSL containers */
//TODO

/* Check our containers */
//TODO
