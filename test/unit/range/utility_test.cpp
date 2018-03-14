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

#include <seqan3/range/utility.hpp>

using namespace seqan3;

struct comparable_range
{
    using reference = std::vector<int>::reference;
    using value_type = std::vector<int>::value_type;

    comparable_range(std::initializer_list<int> list) :
        vec{list}
    {}

    std::vector<int>::iterator begin()
    {
        return vec.begin();
    }

    std::vector<int>::iterator end()
    {
        return vec.end();
    }

    typename std::vector<int>::const_iterator begin() const
    {
        return vec.begin();
    }

    typename std::vector<int>::const_iterator end() const
    {
        return vec.end();
    }

private:
    std::vector<int> vec{};
};

class generic_comparator_operator : public ::testing::Test
{
public:
    comparable_range range_1{{1,2,3,4,5,6}};
    comparable_range range_2{{1,2,3,4,5,6}};
    comparable_range range_3{{1,2,3}};

    comparable_range range_4{{0,1,2,3,4,5}};
};

TEST_F(generic_comparator_operator, equality_comparator)
{
    EXPECT_TRUE(input_range_concept<comparable_range>);

    EXPECT_TRUE(this->range_1 == this->range_2);
    EXPECT_FALSE(this->range_1 == this->range_3);
    EXPECT_FALSE(this->range_2 == this->range_3);
}

TEST_F(generic_comparator_operator, inequality_comparator)
{
    EXPECT_FALSE(this->range_1 != this->range_2);
    EXPECT_TRUE(this->range_1 != this->range_3);
    EXPECT_TRUE(this->range_2 != this->range_3);
}

TEST_F(generic_comparator_operator, less_than_comparator)
{
    EXPECT_FALSE(this->range_1 < this->range_2);
    EXPECT_FALSE(this->range_2 < this->range_1);
    EXPECT_TRUE(this->range_3 < this->range_1);
    EXPECT_TRUE(this->range_4 < this->range_1);
    EXPECT_TRUE(this->range_4 < this->range_3);
}

TEST_F(generic_comparator_operator, less_or_equal_than_comparator)
{
    EXPECT_TRUE(this->range_1 <= this->range_2);
    EXPECT_TRUE(this->range_2 <= this->range_1);
    EXPECT_TRUE(this->range_3 <= this->range_1);
    EXPECT_TRUE(this->range_4 <= this->range_1);
    EXPECT_TRUE(this->range_4 <= this->range_3);
}

TEST_F(generic_comparator_operator, greater_than_comparator)
{
    EXPECT_FALSE(this->range_1 > this->range_2);
    EXPECT_FALSE(this->range_2 > this->range_1);
    EXPECT_FALSE(this->range_3 > this->range_1);
    EXPECT_FALSE(this->range_4 > this->range_1);
    EXPECT_FALSE(this->range_4 > this->range_3);
}

TEST_F(generic_comparator_operator, greater_or_equal_than_comparator)
{
    EXPECT_TRUE(this->range_1 >= this->range_2);
    EXPECT_TRUE(this->range_2 >= this->range_1);
    EXPECT_FALSE(this->range_3 >= this->range_1);
    EXPECT_FALSE(this->range_4 >= this->range_1);
    EXPECT_FALSE(this->range_4 >= this->range_3);
}
