// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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

#include <range/v3/algorithm/equal.hpp>
#include <range/v3/range_traits.hpp>
#include <range/v3/view/take.hpp>
#include <range/v3/view/c_str.hpp>

#include <seqan3/range/container/dummy_container.hpp>
#include <seqan3/range/container/all.hpp>
#include <seqan3/range/view/convert.hpp>

using namespace seqan3;

template <typename T>
class container : public ::testing::Test
{};

TEST(container, concepts)
{
    EXPECT_TRUE(random_access_container_concept<dummy_container<char>>);
}

TEST(container, construction)
{
    dummy_container<char> t1;
    dummy_container<char> t2{};
    EXPECT_EQ(t1.size(), t2.size());

    // initializer list
    dummy_container<char> t3{'A', 'A', 'A', 'A', 'A'};
    dummy_container<char> t4 = {'C', 'C', 'C', 'C', 'C'};
    EXPECT_EQ(t3.size(), t4.size());

    // n * value
    dummy_container<char> t5{2, 'T'};

    // from another dummy_container<char>'s sub-range
    dummy_container<char> t6{t3.begin() + 1, t3.begin() + 3};
    EXPECT_EQ(t5.size(), t6.size());

    // direct from another container
    dummy_container<char> t7{std::string("GGGGG")}; // Note: t7{"GGGGG"} has length 6 because of the \0 character
    EXPECT_EQ(t3.size(), t7.size());
}

TEST(container, assign)
{
    dummy_container<char> t0{'C', 'C'};
    dummy_container<char> t1{'A', 'A', 'A', 'A', 'A'};

    // n * value
    dummy_container<char> t3;
    t3.assign(2, 'C');
    EXPECT_EQ(t3.size(), t0.size());

    // from another container's sub-range
    dummy_container<char> t4;
    t4.assign(t1.cbegin(), t1.cend());
    EXPECT_EQ(t4.size(), t1.size());

    // initializer list
    dummy_container<char> t5, t6;
    t5.assign({'A', 'A', 'A', 'A', 'A'});
    t6 = {'C', 'C', 'C', 'C', 'C'};
    EXPECT_EQ(t5.size(), t1.size());
    EXPECT_EQ(t6.size(), t1.size());

    dummy_container<char> t7;
    t7.assign(std::string("GGGGG")); // Note: t7{"GGGGG"} has length 6 because of the \0 character
    EXPECT_EQ(t7.size(), t1.size());
}

TEST(container, iterators)
{
    dummy_container<char> t1{'A', 'C', 'C', 'G', 'T'};
    dummy_container<char> const t2{'A', 'C', 'C', 'G', 'T'};

    // begin
    EXPECT_THROW(*t1.begin(), std::logic_error);
    EXPECT_THROW(*t1.cbegin(), std::logic_error);
    EXPECT_THROW(*t2.begin(), std::logic_error);
    EXPECT_THROW(*t2.cbegin(), std::logic_error);

    // end and arithmetic
    EXPECT_THROW(*(t1.end()  - 1), std::logic_error);
    EXPECT_THROW(*(t1.cend() - 1), std::logic_error);
    EXPECT_THROW(*(t2.end()  - 1), std::logic_error);
    EXPECT_THROW(*(t2.cend() - 1), std::logic_error);

    // convertibility between const and non-const
    EXPECT_TRUE(t1.cend() == t1.end());

    // range behaviour
    for (auto it = t1.begin(); it != t1.end(); it++)
        EXPECT_THROW(*it, std::logic_error);
}

TEST(container, element_access)
{
    dummy_container<char> t1{'A', 'C', 'C', 'G', 'T'};
    dummy_container<char> const t2{'A', 'C', 'C', 'G', 'T'};

    // at
    EXPECT_THROW(t1.at(0), std::logic_error);
    EXPECT_THROW(t2.at(0), std::logic_error);

    //TODO check at's ability to throw

    // []
    EXPECT_THROW(t1[0], std::logic_error);
    EXPECT_THROW(t2[0], std::logic_error);

    // front
    EXPECT_THROW(t1.front(), std::logic_error);
    EXPECT_THROW(t2.front(), std::logic_error);

    // back
    EXPECT_THROW(t1.back(), std::logic_error);
    EXPECT_THROW(t2.back(), std::logic_error);
}

TEST(container, capacity)
{
    dummy_container<char> t0{};
    dummy_container<char> t1{'A', 'C', 'C', 'G', 'T'};
    dummy_container<char> const t2{'A', 'C', 'C', 'G', 'T'};

    // empty
    EXPECT_TRUE(t0.empty());
    EXPECT_FALSE(t1.empty());
    EXPECT_FALSE(t2.empty());

    // size
    EXPECT_EQ(t0.size(), 0u);
    EXPECT_EQ(t1.size(), 5u);
    EXPECT_EQ(t2.size(), 5u);

    // max_size
    EXPECT_EQ(t0.max_size(), std::numeric_limits<dummy_container<char>::size_type>::max());
    EXPECT_EQ(t1.max_size(), std::numeric_limits<dummy_container<char>::size_type>::max());
    EXPECT_EQ(t2.max_size(), std::numeric_limits<dummy_container<char>::size_type>::max());
}

TEST(container, clear)
{
    dummy_container<char> t0{};
    dummy_container<char> t1{'A', 'C', 'C', 'G', 'T'};

    t1.clear();
    EXPECT_EQ(t0.size(), t1.size());
}

TEST(container, insert)
{
    dummy_container<char> t0{};
    dummy_container<char> t1{'T', 'T', 'T', 'T', 'T'};

    // position, value
    t0.insert(t0.cend(), 'G');
    t0.insert(t0.cend(), 'G');
    t0.insert(t0.cend(), 'G');
    t0.insert(t0.cend(), 'G');
    t0.insert(t0.cbegin() + 1, 'G');
    EXPECT_EQ(t0.size(), t1.size());

    // position, n times values
    t0.clear();
    t0.insert(t0.cend(), 2, 'G');
    t0.insert(t0.cend(), 1, 'G');
    t0.insert(t0.cend(), 1, 'G');
    t0.insert(t0.cbegin(), 1, 'G');
    EXPECT_EQ(t0.size(), t1.size());

    // iterator pair
    t0.clear();
    t0.insert(t0.cend(), t1.begin() + 1, t1.begin() + 3);

    t0.insert(t0.cend(),   t1.cend() - 2, t1.cend());
    t0.insert(t0.cbegin(), t1.cbegin(), t1.cbegin() + 1);
    EXPECT_EQ(t0.size(), t1.size());

    // initializer list
    t0.clear();
    t0.insert(t0.cend(), {'A', 'A', 'A', 'A'});
    t0.insert(t0.cbegin() + 1, 'C');
    EXPECT_EQ(t0.size(), t1.size());
}

TEST(container, erase)
{
    dummy_container<char> t1{'A', 'A', 'A', 'A', 'A'};

    // one element
    t1.erase(t1.begin());
    EXPECT_EQ(t1.size(), (dummy_container<char>{'C', 'C', 'C', 'C'}).size());

    // range
    t1.erase(t1.begin() + 1, t1.begin() + 3);
    EXPECT_EQ(t1.size(), (dummy_container<char>{'G', 'G'}).size());

    // invalid range
    t1.erase(t1.begin() + 1, t1.begin());
    EXPECT_EQ(t1.size(), (dummy_container<char>{'G', 'G'}).size());
}

TEST(container, push_pop)
{
    dummy_container<char> t0{};

    // push_back
    t0.push_back('A');
    EXPECT_EQ(t0.size(),  (dummy_container<char>{'A'}).size());
    t0.push_back('C');
    EXPECT_EQ(t0.size(), (dummy_container<char>{'A', 'C'}).size());

    // pop_back
    t0.pop_back();
    EXPECT_EQ(t0.size(), (dummy_container<char>{'A'}).size());
    t0.pop_back();
    EXPECT_EQ(t0.size(), (dummy_container<char>{}).size());
}

TEST(container, resize)
{
    dummy_container<char> t0{};

    // enlarge without values
    t0.resize(3);
    EXPECT_EQ(t0.size(), (dummy_container<char>{'A', 'A', 'A'}).size());

    // enlarge with value
    t0.resize(5, 'C');
    EXPECT_EQ(t0.size(), (dummy_container<char>{'A', 'A', 'A', 'C', 'C'}).size());

    // shrink with value (no effect)
    t0.resize(4, 'G');
    EXPECT_EQ(t0.size(), (dummy_container<char>{'A', 'A', 'A', 'C'}).size());

    // shrink without value
    t0.resize(2);
    EXPECT_EQ(t0.size(), (dummy_container<char>{'A', 'A'}).size());
}

TEST(container, swap)
{
    dummy_container<char> t0{};
    dummy_container<char> t1{'A', 'C', 'C', 'G', 'T'};

    t0.swap(t1);
    EXPECT_EQ(t0.size(), (dummy_container<char>{'A', 'C', 'C', 'G', 'T'}).size());
    EXPECT_EQ(t1.size(), (dummy_container<char>{}).size());
}
