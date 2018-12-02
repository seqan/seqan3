// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================

#include <sstream>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/alphabet/all.hpp>

using namespace seqan3;
using inner_type = std::vector<dna4>;

class gap_decorator_anchor_set_test_fixture : public ::testing::Test
{
protected:
    std::vector<uint8_t> v_empty{};
    std::vector<uint8_t> const v_const_empty{};
    virtual void SetUp()
    {
        // code here will execute just before the test ensues
    }
};

/*
// concept checks
TEST(gap_decorator_anchor_set, concept_checks)
{
    EXPECT_TRUE((std::RandomAccessIterator<seqan3::detail::random_access_iterator<std::vector<int>>>));
    EXPECT_TRUE((std::RandomAccessIterator<seqan3::detail::random_access_iterator<std::vector<int> const>>));
}

// default constructor
TEST(random_access_iterator_test, default_constructor)
{
    [[maybe_unused]] seqan3::detail::random_access_iterator<std::vector<uint8_t>> it;
    [[maybe_unused]] seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it2;
}

// constructor with empty container reference
TEST(random_access_iterator_test, constructor_ref)
{
    // non-const version
    std::vector<uint8_t> v_empty;
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it(v_empty);
    // const version
    std::vector<uint8_t> const v_const_empty;
    seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it2(v_const_empty);
}

// test constructor call with non-empty container reference and subscript operator
TEST_F(random_access_iterator_test_fixture, constructor_ref2)
{
    // non-const version
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it(v);
    EXPECT_EQ('a', it[0]);
    EXPECT_EQ('t', it[1]);
    // const version
    seqan3::detail::random_access_iterator<std::vector<uint8_t> const> it2(v_const);
    EXPECT_EQ('a', it2[0]);
    EXPECT_EQ('t', it2[1]);
}

// test constructor with non-empty container reference and position offset  and subscript operator
TEST_F(random_access_iterator_test_fixture, constructor_ref3)
{
    // non-const version
    seqan3::detail::random_access_iterator<std::array<long int, 3>> it(a, 1);
    EXPECT_EQ(22, it[0]);
    EXPECT_EQ(33, it[1]);
    // const version
    seqan3::detail::random_access_iterator<std::array<long int, 3> const> it2(a_const, 1);
    EXPECT_EQ(22, it2[0]);
    EXPECT_EQ(33, it2[1]);
}
*/
