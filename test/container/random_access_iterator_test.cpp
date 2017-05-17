// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
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
// Author: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
// ==========================================================================
// Test cases for the random access iterator.
// ==========================================================================

#include <seqan3/container/detail/random_access_iterator.hpp>

#include <gtest/gtest.h>
#include <sstream>
#include <vector>
#include <cassert>

class random_access_iterator_test_fixture: public ::testing::Test {
protected:
    std::vector<uint8_t> *v, *v_au, *v_atz, *v_auvwx, *w, *w_bv;

   random_access_iterator_test_fixture() {}

   virtual void SetUp( ) {
       // code here will execute just before the test ensues
       v = new std::vector<uint8_t>{'a', 't'};
       v_au = new std::vector<uint8_t> {'a', 'u'};
       v_atz = new std::vector<uint8_t> {'a', 't', 'z'};
       v_auvwx = new std::vector<uint8_t> {'a', 'u', 'v', 'w', 'x'};
       w = new std::vector<uint8_t>{'c', 't'};
       w_bv = new std::vector<uint8_t>{'b', 'v'};
   }

   virtual void TearDown( ) {
       delete v, v_au, v_atz, v_auvwx, w, w_bv;
   }

   ~random_access_iterator_test_fixture(){}

};

// default constructor
TEST(random_access_iterator_test, default_constructor)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it;
}

// constructor with empty container reference
TEST(random_access_iterator_test, constructor_ref)
{
    std::vector<uint8_t> v;
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it(v);
}

// constructor with non-empty container reference
TEST_F(random_access_iterator_test_fixture, constructor_ref2)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it(*v);
}

// copy constructor with empty container reference
TEST(random_access_iterator_test, cp_constructor1)
{
    std::vector<uint8_t> v;
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_base(v);
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_derived(it_base);
}

// copy constructor with empty container reference
TEST(random_access_iterator_test, cp_constructor2)
{
    std::vector<uint8_t> v;
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_base(v);
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_derived(it_base);
}

// test assignment construction with empty container reference
TEST(random_access_iterator_test, cp_constructor3)
{
    std::vector<uint8_t> v;
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_base(v);
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_derived = it_base;
}

// explicit desctructor call
TEST(random_access_iterator_test, cp_destructor)
{
    typedef typename seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_type;
    std::vector<uint8_t> v;
    it_type it(v);
    it_type * it_ptr;
    it_ptr->it_type::~it_type();
}

/* equality operators compare internal container positions, not their contents! */
TEST_F(random_access_iterator_test_fixture, equality)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v), it_w(*w);
    EXPECT_TRUE(it_v == it_w);
    ++it_v;
    EXPECT_FALSE(it_v == it_w);
    ++it_w;
    EXPECT_TRUE(it_v == it_w);
}

TEST_F(random_access_iterator_test_fixture, inequality)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v), it_w(*w);
    EXPECT_FALSE(it_v != it_w);
    ++it_v;
    EXPECT_TRUE(it_v != it_w);
}

TEST_F(random_access_iterator_test_fixture, less)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v), it_w(*w);
    EXPECT_FALSE(it_v < it_w);
    ++it_w;
    EXPECT_TRUE(it_v < it_w);
}

TEST_F(random_access_iterator_test_fixture, greater)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v_au), it_w(*w);
    EXPECT_FALSE(it_v > it_w);
    ++it_v;
    EXPECT_TRUE(it_v > it_w);
}

TEST_F(random_access_iterator_test_fixture, leq)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v_atz), it_w(*w);
    EXPECT_TRUE(it_v <= it_w);
    ++it_v;
    EXPECT_FALSE(it_v <= it_w);
}

TEST_F(random_access_iterator_test_fixture, geq)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v_au), it_w(*w);
    EXPECT_TRUE(it_v <= it_w);
    ++it_v;
    EXPECT_FALSE(it_v <= it_w);
}

TEST_F(random_access_iterator_test_fixture, postfix_increment)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v_au), it_w{*w_bv};
    EXPECT_TRUE(it_v == it_w);
    it_v++;
    EXPECT_FALSE(it_v == it_w);
    it_w++;
    EXPECT_TRUE(it_v == it_w);
}

TEST_F(random_access_iterator_test_fixture, prefix_decrement)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v), it_w{*w};
    --++it_v;
    EXPECT_TRUE(it_v == it_w);
}

TEST_F(random_access_iterator_test_fixture, postfix_decrement)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v), it_w{*w};
    (++it_v)--;
    EXPECT_TRUE(it_v == it_w);
}

TEST_F(random_access_iterator_test_fixture, dereference)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v_au);
    EXPECT_TRUE(*it_v == 'a');
    ++it_v;
    EXPECT_TRUE(*it_v == 'u');
}

TEST_F(random_access_iterator_test_fixture, random_access_operator)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v_au);
    EXPECT_TRUE(it_v[0] == 'a');
    EXPECT_TRUE(it_v[1] == 'u');
}

TEST_F(random_access_iterator_test_fixture, operator_plus_assignment)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v_auvwx);
    it_v += 1;
    EXPECT_TRUE(*it_v == 'u');
    it_v += 3;
    EXPECT_TRUE(*it_v == 'x');
}

// test operator+ and self-assignment
TEST_F(random_access_iterator_test_fixture, operator_plus)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v_auvwx);
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_w = it_v + 2;
    EXPECT_TRUE(*it_w == 'v');
    it_v = it_v + 2;
    EXPECT_TRUE(*it_v == 'v');
}

// test operator+ and self-assignment
TEST_F(random_access_iterator_test_fixture, friend_operator_plus)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v_auvwx);
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_w = 2 + it_v;
    EXPECT_TRUE(*it_w == 'v');
    it_v = 2 + it_v;
    EXPECT_TRUE(*it_v == 'v');
}

TEST_F(random_access_iterator_test_fixture, operator_minus_assignment)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v_auvwx);
    it_v += 4;
    it_v -= 3;
    EXPECT_TRUE(*it_v == 'u');
}

TEST_F(random_access_iterator_test_fixture, operator_minus)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v_auvwx);
    it_v = it_v + 4;
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_w = it_v - 2;
    EXPECT_TRUE(*it_w == 'v');
    it_v = it_v - 1;
    EXPECT_TRUE(*it_v == 'w');
}

TEST_F(random_access_iterator_test_fixture, difference)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v_auvwx);
    it_v += 4;
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_w(*w);
    it_w += 1;
    EXPECT_TRUE(3 == (it_v - it_w));
    EXPECT_TRUE(-3 == (it_w - it_v));
}

// test const iterator
TEST(const_random_access_iterator_test, constructors)
{
    const std::vector<uint8_t> v{1,2,3,4};
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>, true> it(v);
}
