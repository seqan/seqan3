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

#include <gtest/gtest.h>

#include <seqan3/alphabet/all.hpp>

using namespace seqan3;

/************** TUPLE INHERITANCE **********************/

template <typename type1, typename type2>
struct test_composition : public cartesian_composition<test_composition<type1, type2>, type1, type2>
{
    using base_t = cartesian_composition<test_composition<type1, type2>, type1, type2>;
    using base_t::base_t;
};

template <typename T>
class cartesian_composition_test: public ::testing::Test
{
public:

    T instance = T{value_1(T{}), value_2(T{})};
    T zero_instance;
    size_t tup_size;

    void SetUp()
    {
        if constexpr(std::is_same_v<T, test_composition<dna4, dna5>>)
        {
           zero_instance = T{dna4{}, dna5{}}; // {A, A}
           tup_size = 2;
        }
        else if(std::is_same_v<T, qualified<dna4, phred42>>)
        {
           zero_instance = T{dna4{}, phred42{}}; // {A, 0}
           tup_size = 2;
        }
    }

    auto value_1(test_composition<dna4, dna5> const &)
    {
        return dna4::G;
    }
    auto value_1(qualified<dna4, phred42> const &)
    {
        return dna4::G;
    }

    auto value_2(test_composition<dna4, dna5> const &)
    {
        return dna5::T;
    }
    auto value_2(qualified<dna4, phred42> const &)
    {
        return phred42{7};
    }

    auto values_to_cmp(test_composition<dna4, dna5> const &)
    {
        return std::make_tuple(/*  v1   v2*/dna4::C, dna5::G,
                               /*==v2  <v2*/dna4::C, dna5::A,
                               /* >v1 ==v2*/dna4::G, dna5::G,
                               /* <v1  >v2*/dna4::A, dna5::T);
    }
    auto values_to_cmp(qualified<dna4, phred42> const &)
    {
        return std::make_tuple(/*  v1   v2*/dna4::C, phred42{6},
                               /*==v2  <v2*/dna4::C, phred42{5},
                               /* >v1 ==v2*/dna4::G, phred42{6},
                               /* <v1  >v2*/dna4::A, phred42{7});
    }
};

using composition_types = ::testing::Types<test_composition<dna4, dna5>,
                                           qualified<dna4, phred42>>;

TYPED_TEST_CASE(cartesian_composition_test, composition_types);

// default/zero construction
TYPED_TEST(cartesian_composition_test, ctr)
{
    [[maybe_unused]] TypeParam t1{};
    EXPECT_EQ(std::tuple_size<TypeParam>::value, TestFixture::tup_size);
}

// initialiser-list initialization
TYPED_TEST(cartesian_composition_test, aggr)
{
    TypeParam t1{};
    TypeParam t2 = TestFixture::instance; // test in fixture to be type independent

    EXPECT_NE(t1, t2);
}

// copy assignment
TYPED_TEST(cartesian_composition_test, cp_assgn)
{
    TypeParam t1 = TestFixture::instance;
    TypeParam t2{};
    TypeParam t3{};

    t2 = t1;
    t3 = t1;
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// zero initialization
TYPED_TEST(cartesian_composition_test, zro)
{
    TypeParam t1 = TestFixture::zero_instance;
    TypeParam t2{};

    EXPECT_EQ(t1, t2);
}

// copy construction
TYPED_TEST(cartesian_composition_test, cp_ctr)
{
    TypeParam t1 = TestFixture::instance;
    TypeParam t2{t1};
    TypeParam t3(t1);

    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move construction
TYPED_TEST(cartesian_composition_test, mv_ctr)
{
    TypeParam t0 = TestFixture::instance;
    TypeParam t1 = TestFixture::instance;
    TypeParam t2{std::move(t1)};

    EXPECT_EQ(t2, t0);

    TypeParam t3(std::move(t2));

    EXPECT_EQ(t3, t0);
}

// move assignment
TYPED_TEST(cartesian_composition_test, mv_assgn)
{
    TypeParam t0 = TestFixture::instance;
    TypeParam t1 = TestFixture::instance;
    TypeParam t2{};
    TypeParam t3{};

    t2 = std::move(t1);

    EXPECT_EQ(t2, t0);

    t3 = std::move(t2);

    EXPECT_EQ(t3, t0);
}

// swap
TYPED_TEST(cartesian_composition_test, swap)
{
    TypeParam t0 = TestFixture::instance;
    TypeParam t1 = TestFixture::instance;
    TypeParam t2{};
    TypeParam t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

// get<0> and get<1>
// get<i> for i > 1 is not tested because of typed_test
TYPED_TEST(cartesian_composition_test, get_i)
{
    TypeParam t0 = TestFixture::instance;

    static_assert(std::is_same_v<decltype(seqan3::get<0>(t0)), decltype(TestFixture::value_1(TypeParam{})) &>);
    static_assert(std::is_same_v<decltype(seqan3::get<1>(t0)), decltype(TestFixture::value_2(TypeParam{})) &>);

    EXPECT_EQ(seqan3::get<0>(t0), TestFixture::value_1(TypeParam{}));
    EXPECT_EQ(seqan3::get<1>(t0), TestFixture::value_2(TypeParam{}));
}

// std::get<1>
TYPED_TEST(cartesian_composition_test, stdget_i)
{
    TypeParam t0 = TestFixture::instance;

    static_assert(std::is_same_v<decltype(std::get<0>(t0)), decltype(TestFixture::value_1(TypeParam{})) &>);
    static_assert(std::is_same_v<decltype(std::get<1>(t0)), decltype(TestFixture::value_2(TypeParam{})) &>);

    EXPECT_EQ(std::get<0>(t0), TestFixture::value_1(TypeParam{}));
    EXPECT_EQ(std::get<1>(t0), TestFixture::value_2(TypeParam{}));
}

// structured bindings
TYPED_TEST(cartesian_composition_test, struct_binding)
{
  TypeParam t0 = TestFixture::instance;
  auto [ i, l ] = t0;

  static_assert(std::is_same_v<decltype(i), decltype(TestFixture::value_1(TypeParam{}))>);
  static_assert(std::is_same_v<decltype(l), decltype(TestFixture::value_2(TypeParam{}))>);

  EXPECT_EQ(i, TestFixture::value_1(TypeParam{}));
  EXPECT_EQ(l, TestFixture::value_2(TypeParam{}));
}

// get<type>
TYPED_TEST(cartesian_composition_test, get_type)
{
    TypeParam t0 = TestFixture::instance;

    EXPECT_EQ(seqan3::get<decltype(TestFixture::value_1(TypeParam{}))>(t0), TestFixture::value_1(TypeParam{}));
    EXPECT_EQ(seqan3::get<decltype(TestFixture::value_2(TypeParam{}))>(t0), TestFixture::value_2(TypeParam{}));
}

// std::get<type>
TYPED_TEST(cartesian_composition_test, stdget_type)
{
    TypeParam t0 = TestFixture::instance;

    EXPECT_EQ(std::get<decltype(TestFixture::value_1(TypeParam{}))>(t0), TestFixture::value_1(TypeParam{}));
    EXPECT_EQ(std::get<decltype(TestFixture::value_2(TypeParam{}))>(t0), TestFixture::value_2(TypeParam{}));
}

// Custom constructor that assigns one type and defaults the other values
// (after get<> was tested)
TYPED_TEST(cartesian_composition_test, custom_ctr)
{
    // first type
    TypeParam t1{TestFixture::value_1(TypeParam{})};
    TypeParam t2 = TestFixture::zero_instance;

    EXPECT_NE(get<0>(t1), get<0>(t2));
    EXPECT_EQ(get<1>(t1), get<1>(t2));
    EXPECT_EQ(get<0>(t1), TestFixture::value_1(TypeParam{}));

    // second type
    TypeParam t3{TestFixture::value_2(TypeParam{})};

    EXPECT_EQ(get<0>(t3), get<0>(t2));
    EXPECT_NE(get<1>(t3), get<1>(t2));
    EXPECT_EQ(get<1>(t3), TestFixture::value_2(TypeParam{}));

}

// std::tuple_element
TYPED_TEST(cartesian_composition_test, tuple_element)
{
    using pt = TypeParam;

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, decltype(TestFixture::value_1(TypeParam{}))>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, decltype(TestFixture::value_2(TypeParam{}))>);
}

// type deduction
TYPED_TEST(cartesian_composition_test, type_deduce)
{
    TypeParam t0 = TestFixture::instance;
    using pt = decltype(t0);

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, decltype(TestFixture::value_1(TypeParam{}))>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, decltype(TestFixture::value_2(TypeParam{}))>);
}

// explicit cast to element
TYPED_TEST(cartesian_composition_test, cast_to_element)
{
    TypeParam t0 = TestFixture::instance;

    auto d = static_cast<decltype(TestFixture::value_1(TypeParam{}))>(t0);
    auto q = static_cast<decltype(TestFixture::value_2(TypeParam{}))>(t0);
    static_assert(std::is_same_v<decltype(d), decltype(TestFixture::value_1(TypeParam{}))>);
    static_assert(std::is_same_v<decltype(q), decltype(TestFixture::value_2(TypeParam{}))>);

    EXPECT_EQ(d, TestFixture::value_1(TypeParam{}));
    EXPECT_EQ(q, TestFixture::value_2(TypeParam{}));
}

// comparison operators
TYPED_TEST(cartesian_composition_test, cmp)
{
    TypeParam t0 = {std::get<0>(TestFixture::values_to_cmp(TypeParam{})),
                    std::get<1>(TestFixture::values_to_cmp(TypeParam{}))};
    TypeParam t1 = {std::get<2>(TestFixture::values_to_cmp(TypeParam{})),
                    std::get<3>(TestFixture::values_to_cmp(TypeParam{}))};
    TypeParam t2 = {std::get<4>(TestFixture::values_to_cmp(TypeParam{})),
                    std::get<5>(TestFixture::values_to_cmp(TypeParam{}))};
    TypeParam t3 = {std::get<6>(TestFixture::values_to_cmp(TypeParam{})),
                    std::get<7>(TestFixture::values_to_cmp(TypeParam{}))};

    EXPECT_EQ(t1, t1);

    EXPECT_NE(t0, t1);
    EXPECT_NE(t0, t2);
    EXPECT_NE(t2, t3);

    EXPECT_LT(t0, t2);
    EXPECT_LT(t1, t0);
    EXPECT_LT(t1, t2);
    EXPECT_LT(t3, t0);
    EXPECT_LT(t3, t1);
    EXPECT_LT(t3, t2);

    EXPECT_LE(t0, t2);
    EXPECT_LE(t1, t0);
    EXPECT_LE(t1, t2);
    EXPECT_LE(t3, t0);
    EXPECT_LE(t3, t1);
    EXPECT_LE(t3, t2);
    EXPECT_LE(t1, t1);

    EXPECT_GE(t0, t1);
    EXPECT_GE(t0, t3);
    EXPECT_GE(t1, t1);
    EXPECT_GE(t1, t3);
    EXPECT_GE(t2, t0);
    EXPECT_GE(t2, t1);
    EXPECT_GE(t2, t3);

    EXPECT_GT(t0, t1);
    EXPECT_GT(t0, t3);
    EXPECT_GT(t1, t3);
    EXPECT_GT(t2, t0);
    EXPECT_GT(t2, t1);
    EXPECT_GT(t2, t3);
}
