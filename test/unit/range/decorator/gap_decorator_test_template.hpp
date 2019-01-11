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
//#include <seqan3/range/decorator/all.hpp>
#include <seqan3/range/decorator/gap_decorator_anchor_set.hpp>

using namespace seqan3;
using inner_type = std::vector<dna4>;

template <typename T>
class gap_decorator : public ::testing::Test
{};

// TODO: add more container types
//typedef ::testing::Types<std::vector<dna4>, std::vector<dna15>> ContainerTypes;

TYPED_TEST_CASE_P(gap_decorator);

// concept test
TYPED_TEST_P(gap_decorator, concepts)
{
    EXPECT_TRUE(aligned_sequence_concept<TypeParam>);
}

// default constructor
TYPED_TEST_P(gap_decorator, default_constructor)
{
    [[maybe_unused]] TypeParam  gd;
    //[[maybe_unused]] TypeParam gd_const;
}

/*
// copy constructor
TYPED_TEST_P(gap_decorator, copy_constructor)
{
    TypeParam gd;
    [[maybe_unused]] TypeParam gd_cp(gd);
}

// copy construction via assignment
TYPED_TEST_P(gap_decorator, copy_assign_constructor)
{
    TypeParam gd;
    [[maybe_unused]] TypeParam gd_cp = gd;
}

// move constructor
TYPED_TEST_P(gap_decorator, move_constructor)
{
    TypeParam gd;
    [[maybe_unused]] TypeParam gd_cp(std::move(gd));
    //EXPECT_EQ(gd_cp, gd);
}

// move assignment
TYPED_TEST_P(gap_decorator, move_assignment)
{
    TypeParam gd;
    [[maybe_unused]] TypeParam gd_cp = std::move(gd);
    //EXPECT_EQ(gd_cp, gd);
}

// seqan3::aligned_sequence_concept test


//!\brief Use default deconstructor.
~gap_decorator_anchor_set() = default;
//!\brief Construct by host and explicit position.
constexpr gap_decorator_anchor_set(inner_type * sequence): data{new data_t{sequence}} {};
//!\}
*/
/*
TYPED_TEST_P(alphabet, default_value_constructor)
{
    [[maybe_unused]] TypeParam t1;
    [[maybe_unused]] TypeParam t2{};
}
*/
/*
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

REGISTER_TYPED_TEST_CASE_P(gap_decorator, concepts, default_constructor);
//, copy_constructor, copy_assign_constructor, move_constructor, move_assignment);
