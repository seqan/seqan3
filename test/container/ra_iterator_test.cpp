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

#include <seqan3/container/detail/ra_iterator.hpp>

#include <gtest/gtest.h>
#include <sstream>
#include <vector>

//#include <seqan3/container/concepts.hpp>
//#include <seqan3/alphabet/alphabet_container.hpp>

/* test container concept properties */
// default constructor
TEST(random_access_iterator_test, default_constructor)
{
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it;
}

// constructor with empty container reference
TEST(random_access_iterator_test, constructor_ref)
{
    std::vector<uint8_t> v;
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it(v);
}

// constructor with non-empty container reference
TEST(random_access_iterator_test, constructor_ref2)
{
    std::vector<uint8_t> v{'a', 't'};
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it(v);
}

// copy constructor with empty container reference
TEST(random_access_iterator_test, cp_constructor1)
{
    std::vector<uint8_t> v;
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_base(v);
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_derived(it_base);
}

// copy constructor with empty container reference
TEST(random_access_iterator_test, cp_constructor2)
{
    std::vector<uint8_t> v;
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_base(v);
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_derived(it_base);
}

// TODO: test for mv constructor ra_iterator(ra_iterator &&) = default;

// test assignment construction with empty container reference
TEST(random_access_iterator_test, cp_constructor3)
{
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_base;
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_derived = it_base;
}

// explicit desctructor call
TEST(random_access_iterator_test, cp_destructor)
{
    typedef typename seqan3::detail::ra_iterator<std::vector<uint8_t>> it_type;
    std::vector<uint8_t> v;
    it_type it(v);
    it_type * it_ptr;
    it_ptr->it_type::~it_type();
}

/* equality operators */
TEST(random_access_iterator_test, equality)
{
    std::vector<uint8_t> v{'a', 't'}, w{'c', 't'};
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_v(v), it_w(w);
    EXPECT_FALSE(*it_v == *it_w);
    EXPECT_TRUE(it_v[1] == it_w[1]);
}

TEST(random_access_iterator_test, inequality)
{
    std::vector<uint8_t> v{'a', 't'}, w{'c', 't'};
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_v(v), it_w(w);
    EXPECT_TRUE(*it_v != *it_w);
    EXPECT_FALSE(it_v[1] != it_w[1]);
}

TEST(random_access_iterator_test, less)
{
    std::vector<uint8_t> v{'a', 't'}, w{'c', 't'};
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_v(v), it_w(w);
    EXPECT_TRUE(*it_v < *it_w);
    EXPECT_FALSE(it_v[1] < it_w[1]);
}

TEST(random_access_iterator_test, greater)
{
    std::vector<uint8_t> v{'a', 'u'}, w{'c', 't'};
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_v(v), it_w(w);
    EXPECT_FALSE(*it_v > *it_w);
    EXPECT_TRUE(it_v[1] > it_w[1]);
}

TEST(random_access_iterator_test, leq)
{
    std::vector<uint8_t> v{'a', 't', 'z'}, w{'c', 't'};
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_v(v), it_w(w);
    EXPECT_TRUE(*it_v <= *it_w);
    EXPECT_TRUE(it_v[1] <= it_w[1]);
    EXPECT_FALSE(it_v[2] <= it_w[1]);
}

TEST(random_access_iterator_test, geq)
{
    std::vector<uint8_t> v{'a', 'u'}, w{'c', 't'};
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_v(v), it_w(w);
    EXPECT_TRUE(*it_v <= *it_w);
    EXPECT_FALSE(it_v[1] <= it_w[1]);
}

TEST(random_access_iterator_test, prefix_increment)
{
    std::vector<uint8_t> v{'a', 'u'};
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_v(v);
    ++it_v;
    EXPECT_EQ(*it_v, 'u');
}

TEST(random_access_iterator_test, postfix_increment)
{
    std::vector<uint8_t> v{'a', 'u'};
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_v(v);
    it_v++;
    EXPECT_EQ(*it_v, 'u');
}

TEST(random_access_iterator_test, prefix_decrement)
{
    std::vector<uint8_t> v{'a', 'u'};
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_v(v);
    --++it_v;
    EXPECT_EQ(*it_v, 'a');
}

TEST(random_access_iterator_test, postfix_decrement)
{
    std::vector<uint8_t> v{'a', 'u'};
    seqan3::detail::ra_iterator<std::vector<uint8_t>> it_v(v);
    (++it_v)--;
    EXPECT_EQ(*it_v, 'a');
}
