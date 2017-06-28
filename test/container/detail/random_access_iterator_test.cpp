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
// ============================================================================

/*!\file container/detail/random_access_iterator_test.hpp
 * \brief Random access iterator for const and non const containers.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \ingroup iterator test
 */

#include <seqan3/container/detail/random_access_iterator.hpp>

#include <gtest/gtest.h>
#include <sstream>
#include <vector>

class random_access_iterator_test_fixture: public ::testing::Test {
protected:

   std::vector<uint8_t> v_empty;
   const std::vector<uint8_t> v_const_empty{};
   std::vector<uint8_t> *v, *v2, *v3, *v4, *w, *w2;
   const std::vector<uint8_t> *v_const, *v2_const, *v3_const, *v4_const, *w_const, *w2_const;
   std::array<long int, 3> a, b;
   const std::array<long int, 3> a_const[3] = {11, 22, 33}; //, b_const;

   virtual void SetUp( ) {
       // code here will execute just before the test ensues
       v = new std::vector<uint8_t>                 {'a', 't'};
       v_const = new const std::vector<uint8_t>     {'a', 't'};
       v2 = new std::vector<uint8_t>                {'a', 'u'};
       v2_const = new const std::vector<uint8_t>    {'a', 'u'};
       v3 = new std::vector<uint8_t>                {'a', 't', 'z'};
       v3_const = new const std::vector<uint8_t>    {'a', 't', 'z'};
       v4 = new std::vector<uint8_t>                {'a', 'u', 'v', 'w', 'x'};
       v4_const = new const std::vector<uint8_t>    {'a', 'u', 'v', 'w', 'x'};
       w = new std::vector<uint8_t>                 {'c', 't'};
       w_const = new const std::vector<uint8_t>     {'c', 't'};
       w2 = new std::vector<uint8_t>                {'b', 'v'};
       w2_const = new const std::vector<uint8_t>    {'b', 'v'};
       a = {11, 22, 33};
       //a_const = {11, 22, 33};
   }

   virtual void TearDown( ) {
       delete v, delete v2, delete v3, delete v4, delete w, delete w2;
       delete v_const, delete v2_const, delete v3_const, delete v4_const, delete w_const, delete w2_const;
   }

   ~random_access_iterator_test_fixture(){}
};

// default constructor
TEST(random_access_iterator_test, default_constructor)
{
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it;
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it2;
}

// constructor with empty container reference
TEST(random_access_iterator_test, constructor_ref)
{
    // non-const version
    std::vector<uint8_t> v_empty;
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it(v_empty);
    // const version
    const std::vector<uint8_t> v_const_empty;
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it2(v_const_empty);
}

// test constructor call with non-empty container reference and subscript operator
TEST_F(random_access_iterator_test_fixture, constructor_ref2)
{
    // non-const version
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it(*v);
    EXPECT_EQ('a', it[0]);
    EXPECT_EQ('t', it[1]);
    // const version
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it2(*v_const);
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
    seqan3::detail::random_access_iterator<const std::array<long int, 3>> it2(*a_const, 1);
    EXPECT_EQ(22, it2[0]);
    EXPECT_EQ(33, it2[1]);
}

// copy constructor with empty container reference
TEST_F(random_access_iterator_test_fixture, cp_constructor1)
{
    // non-const container
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_base(v_empty);
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_derivate(it_base);
    // const container
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_base2(v_const_empty);
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_derivate2(it_base2);
}

// copy constructor with non-empty container reference
TEST_F(random_access_iterator_test_fixture, cp_constructor2)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_base(*v);
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_derivate(it_base);
    EXPECT_EQ('a', it_base[0]);
    EXPECT_EQ('a', it_derivate[0]);
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_base2(*v_const);
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_derivate2(it_base2);
    EXPECT_EQ('a', it_base2[0]);
    EXPECT_EQ('a', it_derivate2[0]);

}

// test assignment construction with empty container reference
TEST_F(random_access_iterator_test_fixture, constructor_assign1)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_base(v_empty);
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_derived = it_base;
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_base2(v_const_empty);
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_derived2 = it_base2;
}

// test assignment construction with non-empty container reference and subscript
TEST_F(random_access_iterator_test_fixture, constructor_assign2)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_base(*v);
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_derivate = it_base;
    EXPECT_EQ('t', it_base[1]);
    EXPECT_EQ('t', it_derivate[1]);
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_base2(*v_const);
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_derivate2 = it_base2;
    EXPECT_EQ('t', it_base2[1]);
    EXPECT_EQ('t', it_derivate2[1]);
}

// test move constructor
TEST_F(random_access_iterator_test_fixture, constructor_move)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it1(*v);
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it2(std::move(it1));
    EXPECT_EQ('a', it2[0]);
    EXPECT_EQ('t', it2[1]);
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it3(*v_const);
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it4(std::move(it3));
    EXPECT_EQ('a', it3[0]);
    EXPECT_EQ('t', it4[1]);
}

// test move assignment
TEST_F(random_access_iterator_test_fixture, move_assign)
{
    // non-const
    seqan3::detail::random_access_iterator<std::array<long int, 3>> it1, it2;
    it2 = std::move(it1);
    // const
    seqan3::detail::random_access_iterator<const std::array<long int, 3>> it3, it4;
    it4 = std::move(it3);
}

// explicit desctructor call
TEST_F(random_access_iterator_test_fixture, cp_destructor)
{
    // non-const
    using iterator_type = typename seqan3::detail::random_access_iterator<std::vector<uint8_t>>;
    iterator_type it(v_empty);
    iterator_type * it_ptr;
    it_ptr->iterator_type::~iterator_type();
    // const
    using iterator_type2 = typename seqan3::detail::random_access_iterator<const std::vector<uint8_t>>;
    iterator_type2 it2(v_const_empty);
    iterator_type2 * it_ptr2;
    it_ptr2->iterator_type2::~iterator_type2();

}

// equality operators compare internal container positions, not their contents!
TEST_F(random_access_iterator_test_fixture, equality)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v), it_w(*w);
    EXPECT_TRUE(it_v == it_w);
    ++it_v;
    EXPECT_FALSE(it_v == it_w);
    ++it_w;
    EXPECT_TRUE(it_v == it_w);
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v_const), it_w2(*w_const);
    EXPECT_TRUE(it_v2 == it_w2);
    ++it_v2;
    EXPECT_FALSE(it_v2 == it_w2);
    ++it_w2;
    EXPECT_TRUE(it_v2 == it_w2);
}

TEST_F(random_access_iterator_test_fixture, inequality)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v), it_w(*w);
    EXPECT_FALSE(it_v != it_w);
    ++it_v;
    EXPECT_TRUE(it_v != it_w);
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v_const), it_w2(*w_const);
    EXPECT_FALSE(it_v2 != it_w2);
    ++it_v2;
    EXPECT_TRUE(it_v2 != it_w2);
}

TEST_F(random_access_iterator_test_fixture, less)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v), it_w(*w);
    EXPECT_FALSE(it_v < it_w);
    ++it_w;
    EXPECT_TRUE(it_v < it_w);
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v_const), it_w2(*w_const);
    EXPECT_FALSE(it_v2 < it_w2);
    ++it_w2;
    EXPECT_TRUE(it_v2 < it_w2);
}

TEST_F(random_access_iterator_test_fixture, greater)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v2), it_w(*w);
    EXPECT_FALSE(it_v > it_w);
    ++it_v;
    EXPECT_TRUE(it_v > it_w);
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v2_const), it_w2(*w_const);
    EXPECT_FALSE(it_v2 > it_w2);
    ++it_v2;
    EXPECT_TRUE(it_v2 > it_w2);

}

TEST_F(random_access_iterator_test_fixture, leq)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v3), it_w(*w);
    EXPECT_TRUE(it_v <= it_w);
    ++it_v;
    EXPECT_FALSE(it_v <= it_w);
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v3_const), it_w2(*w_const);
    EXPECT_TRUE(it_v2 <= it_w2);
    ++it_v2;
    EXPECT_FALSE(it_v2 <= it_w2);
}


TEST_F(random_access_iterator_test_fixture, geq)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v2), it_w(*w);
    EXPECT_TRUE(it_v <= it_w);
    ++it_v;
    EXPECT_FALSE(it_v <= it_w);
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v2_const), it_w2(*w_const);
    EXPECT_TRUE(it_v2 <= it_w2);
    ++it_v2;
    EXPECT_FALSE(it_v2 <= it_w2);
}

TEST_F(random_access_iterator_test_fixture, postfix_increment)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v2), it_w{*w2};
    EXPECT_TRUE(it_v == it_w);
    it_v++;
    EXPECT_FALSE(it_v == it_w);
    it_w++;
    EXPECT_TRUE(it_v == it_w);
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v2_const), it_w2{*w2_const};
    EXPECT_TRUE(it_v2 == it_w2);
    it_v2++;
    EXPECT_FALSE(it_v2 == it_w2);
    it_w2++;
    EXPECT_TRUE(it_v2 == it_w2);
}

TEST_F(random_access_iterator_test_fixture, prefix_decrement)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v), it_w{*w};
    --++it_v;
    EXPECT_TRUE(it_v == it_w);
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v_const), it_w2{*w_const};
    --++it_v2;
    EXPECT_TRUE(it_v2 == it_w2);
}


TEST_F(random_access_iterator_test_fixture, postfix_decrement)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v), it_w{*w};
    (++it_v)--;
    EXPECT_TRUE(it_v == it_w);
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v_const), it_w2{*w_const};
    (++it_v2)--;
    EXPECT_TRUE(it_v2 == it_w2);
}

TEST_F(random_access_iterator_test_fixture, dereference)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v2);
    EXPECT_TRUE(*it_v == 'a');
    ++it_v;
    EXPECT_TRUE(*it_v == 'u');
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v2_const);
    EXPECT_TRUE(*it_v2 == 'a');
    ++it_v2;
    EXPECT_TRUE(*it_v2 == 'u');
}

TEST_F(random_access_iterator_test_fixture, random_access_operator)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v2);
    EXPECT_TRUE(it_v[0] == 'a');
    EXPECT_TRUE(it_v[1] == 'u');
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v2_const);
    EXPECT_TRUE(it_v2[0] == 'a');
    EXPECT_TRUE(it_v2[1] == 'u');
}

TEST_F(random_access_iterator_test_fixture, operator_plus_assignment)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v4);
    it_v += 1;
    EXPECT_TRUE(*it_v == 'u');
    it_v += 3;
    EXPECT_TRUE(*it_v == 'x');
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v4_const);
    it_v2 += 1;
    EXPECT_TRUE(*it_v2 == 'u');
    it_v2 += 3;
    EXPECT_TRUE(*it_v2 == 'x');
}

// test operator+ and self-assignment
TEST_F(random_access_iterator_test_fixture, operator_plus)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v4);
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_w = it_v + 2;
    EXPECT_TRUE(*it_w == 'v');
    it_v = it_v + 2;
    EXPECT_TRUE(*it_v == 'v');
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v4_const);
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_w2 = it_v2 + 2;
    EXPECT_TRUE(*it_w2 == 'v');
    it_v2 = it_v2 + 2;
    EXPECT_TRUE(*it_v2 == 'v');
}

// test operator+ and self-assignment
TEST_F(random_access_iterator_test_fixture, friend_operator_plus)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v4);
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_w = 2 + it_v;
    EXPECT_TRUE(*it_w == 'v');
    it_v = 2 + it_v;
    EXPECT_TRUE(*it_v == 'v');
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v4_const);
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_w2 = 2 + it_v2;
    EXPECT_TRUE(*it_w2 == 'v');
    it_v2 = 2 + it_v2;
    EXPECT_TRUE(*it_v2 == 'v');
}

TEST_F(random_access_iterator_test_fixture, operator_minus_assignment)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v4);
    it_v += 4;
    it_v -= 3;
    EXPECT_TRUE(*it_v == 'u');
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v4_const);
    it_v2 += 4;
    it_v2 -= 3;
    EXPECT_TRUE(*it_v2 == 'u');
}

TEST_F(random_access_iterator_test_fixture, operator_minus)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v4);
    it_v = it_v + 4;
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_w = it_v - 2;
    EXPECT_TRUE(*it_w == 'v');
    it_v = it_v - 1;
    EXPECT_TRUE(*it_v == 'w');
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v4_const);
    it_v2 = it_v2 + 4;
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_w2 = it_v2 - 2;
    EXPECT_TRUE(*it_w2 == 'v');
    it_v2 = it_v2 - 1;
    EXPECT_TRUE(*it_v2 == 'w');
}

TEST_F(random_access_iterator_test_fixture, difference)
{
    // non-const
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_v(*v4);
    it_v += 4;
    seqan3::detail::random_access_iterator<std::vector<uint8_t>> it_w(*w);
    it_w += 1;
    EXPECT_TRUE(3 == (it_v - it_w));
    EXPECT_TRUE(-3 == (it_w - it_v));
    // const
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_v2(*v4_const);
    it_v2 += 4;
    seqan3::detail::random_access_iterator<const std::vector<uint8_t>> it_w2(*w_const);
    it_w2 += 1;
    EXPECT_TRUE(3 == (it_v2 - it_w2));
    EXPECT_TRUE(-3 == (it_w2 - it_v2));
}
