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

template <typename type1, typename type2>
struct test_composition : public cartesian_composition<test_composition<type1, type2>, type1, type2>
{
    using base_t = cartesian_composition<test_composition<type1, type2>, type1, type2>;
    using base_t::base_t;
    using base_t::operator=;

    using base_t::operator==;
    using base_t::operator!=;
    using base_t::operator<;
    using base_t::operator>;
    using base_t::operator<=;
    using base_t::operator>=;
};

template <typename T>
class cartesian_composition_test : public ::testing::Test {};

template <>
class cartesian_composition_test<test_composition<dna4, dna5>> : public ::testing::Test
{
public:
    using T = test_composition<dna4, dna5>;

    T instance = T{value_1(), value_2()};
    T zero_instance = T{decltype(value_1()){}, decltype(value_2()){}};
    size_t tup_size{2};

    // test_composition<dna4, dna5>
    // -------------------------------------------------------------------------
    dna4 value_1()
    {
        return dna4::G;
    }
    rna4 assignable_to_value_1()
    {
        return rna4::G;
    }
    dna5 value_2()
    {
        return dna5::G;
    }
    rna5 assignable_to_value_2()
    {
        return rna5::G;
    }
    auto values_to_cmp()
    {
        return std::make_tuple(/*low */dna4::A, dna5::A,
                               /*mid */dna4::C, dna5::C,
                               /*high*/dna4::T, dna5::T);
    }
};

template <>
class cartesian_composition_test<qualified<dna4, phred42>>: public ::testing::Test
{
public:
    using T = qualified<dna4, phred42>;

    T instance = T{value_1(), value_2()};
    T zero_instance = T{decltype(value_1()){}, decltype(value_2()){}};
    size_t tup_size{2};

    // qualified<dna4, phred42>
    // -------------------------------------------------------------------------
    dna4 value_1()
    {
        return dna4::G;
    }
    rna4 assignable_to_value_1()
    {
        return rna4::G;
    }
    phred42 value_2()
    {
        return phred42{6};
    }
    phred42 assignable_to_value_2()
    {
        return phred42{6}; // replace if assignable subtype becomes available
    }
    auto values_to_cmp()
    {
        return std::make_tuple(/*low */dna4::A, phred42{1},
                               /*mid */dna4::C, phred42{4},
                               /*high*/dna4::T, phred42{9});
    }
};

template <>
class cartesian_composition_test<structured_rna<rna4, dot_bracket3>>: public ::testing::Test
{
public:
    using T = structured_rna<rna4, dot_bracket3>;

    T instance = T{value_1(), value_2()};
    T zero_instance = T{decltype(value_1()){}, decltype(value_2()){}};
    size_t tup_size{2};

    // structured_rna<rna4, dot_bracket3>
    // -------------------------------------------------------------------------
    rna4 value_1()
    {
        return rna4::G;
    }
    dna4 assignable_to_value_1()
    {
        return dna4::G;
    }
    dot_bracket3 value_2()
    {
        return dot_bracket3::PAIR_OPEN;
    }
    dot_bracket3 assignable_to_value_2()
    {
        return dot_bracket3::PAIR_OPEN; // replace if assignable subtype becomes available
    }
    auto values_to_cmp()
    {
        return std::make_tuple(/*low */rna4::A, dot_bracket3::UNPAIRED,
                               /*mid */rna4::C, dot_bracket3::PAIR_OPEN,
                               /*high*/rna4::T, dot_bracket3::PAIR_CLOSE);
    }

};

template <>
class cartesian_composition_test<structured_aa<aa27, dssp9>>: public ::testing::Test
{
public:
    using T = structured_aa<aa27, dssp9>;

    T instance = T{value_1(), value_2()};
    T zero_instance = T{decltype(value_1()){}, decltype(value_2()){}};
    size_t tup_size{2};

    // structured_aa<aa27, dssp9>
    // -------------------------------------------------------------------------
    aa27 value_1()
    {
        return aa27::K;
    }
    aa27 assignable_to_value_1()
    {
        return aa27::K; // replace if assignable subtype becomes available
    }
    dssp9 value_2()
    {
        return dssp9::I;
    }
    dssp9 assignable_to_value_2()
    {
        return dssp9::I; // replace if assignable subtype becomes available
    }
    auto values_to_cmp()
    {
        return std::make_tuple(/*low */aa27::A, dssp9::H,
                               /*mid */aa27::P, dssp9::I,
                               /*high*/aa27::Z, dssp9::X);
    }
};

template <>
class cartesian_composition_test<cigar<>>: public ::testing::Test
{
public:
    using T = cigar<>;

    T instance = T{value_1(), value_2()};
    T zero_instance = T{decltype(value_1()){}, decltype(value_2()){}};
    size_t tup_size{2};

    // cigar
    // -------------------------------------------------------------------------
    cigar_op value_1()
    {
        return cigar_op::D;
    }
    cigar_op assignable_to_value_1()
    {
        return cigar_op::D; // replace if assignable subtype becomes available
    }
    uint32_t value_2()
    {
        return 200u;
    }
    uint32_t assignable_to_value_2()
    {
        return uint8_t(200); // replace if assignable subtype becomes available
    }
    auto values_to_cmp()
    {
        return std::make_tuple(/*low */cigar_op::M,  (uint32_t)1u,
                               /*mid */cigar_op::X,  (uint32_t)100u,
                               /*high*/cigar_op::EQ, (uint32_t)1000u);
    }
};

using composition_types = ::testing::Types<test_composition<dna4, dna5>,
                                           structured_rna<rna4, dot_bracket3>,
                                           structured_aa<aa27, dssp9>,
                                           qualified<dna4, phred42>,
                                           cigar<>>;

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

    static_assert(std::is_same_v<decltype(seqan3::get<0>(t0)), decltype(TestFixture::value_1()) &>);
    static_assert(std::is_same_v<decltype(seqan3::get<1>(t0)), decltype(TestFixture::value_2()) &>);

    EXPECT_EQ(seqan3::get<0>(t0), TestFixture::value_1());
    EXPECT_EQ(seqan3::get<1>(t0), TestFixture::value_2());
}

// std::get<1>
TYPED_TEST(cartesian_composition_test, stdget_i)
{
    TypeParam t0 = TestFixture::instance;

    static_assert(std::is_same_v<decltype(std::get<0>(t0)), decltype(TestFixture::value_1()) &>);
    static_assert(std::is_same_v<decltype(std::get<1>(t0)), decltype(TestFixture::value_2()) &>);

    EXPECT_EQ(std::get<0>(t0), TestFixture::value_1());
    EXPECT_EQ(std::get<1>(t0), TestFixture::value_2());
}

// structured bindings
TYPED_TEST(cartesian_composition_test, struct_binding)
{
  TypeParam t0 = TestFixture::instance;
  auto [ i, l ] = t0;

  static_assert(std::is_same_v<decltype(i), decltype(TestFixture::value_1())>);
  static_assert(std::is_same_v<decltype(l), decltype(TestFixture::value_2())>);

  EXPECT_EQ(i, TestFixture::value_1());
  EXPECT_EQ(l, TestFixture::value_2());
}

// get<type>
TYPED_TEST(cartesian_composition_test, get_type)
{
    TypeParam t0 = TestFixture::instance;

    EXPECT_EQ(seqan3::get<decltype(TestFixture::value_1())>(t0), TestFixture::value_1());
    EXPECT_EQ(seqan3::get<decltype(TestFixture::value_2())>(t0), TestFixture::value_2());
}

// std::get<type>
TYPED_TEST(cartesian_composition_test, stdget_type)
{
    TypeParam t0 = TestFixture::instance;

    EXPECT_EQ(std::get<decltype(TestFixture::value_1())>(t0), TestFixture::value_1());
    EXPECT_EQ(std::get<decltype(TestFixture::value_2())>(t0), TestFixture::value_2());
}

// Custom constructor that assigns one type and defaults the other values
// (after get<> was tested)
TYPED_TEST(cartesian_composition_test, custom_ctr)
{
    // first type
    TypeParam t1{TestFixture::value_1()};
    TypeParam t2 = TestFixture::zero_instance;

    EXPECT_NE(get<0>(t1), get<0>(t2));
    EXPECT_EQ(get<1>(t1), get<1>(t2));
    EXPECT_EQ(get<0>(t1), TestFixture::value_1());

    // second type
    TypeParam t3{TestFixture::value_2()};

    EXPECT_EQ(get<0>(t3), get<0>(t2));
    EXPECT_NE(get<1>(t3), get<1>(t2));
    EXPECT_EQ(get<1>(t3), TestFixture::value_2());

}

// Custom constructor that assigns one type from an assignable subtype and defaults the other values
TYPED_TEST(cartesian_composition_test, custom_ctr_subtype)
{
    // first type
    TypeParam t1{TestFixture::assignable_to_value_1()};
    TypeParam t_d{};

    EXPECT_EQ(get<0>(t1), TestFixture::value_1());
    EXPECT_EQ(get<1>(t1), get<1>(t_d));

    // second type
    TypeParam t3{TestFixture::assignable_to_value_2()};

    EXPECT_EQ(get<0>(t3), get<0>(t_d));
    EXPECT_EQ(get<1>(t3), TestFixture::value_2());
}

// Custom assignment operator that assigns one type and defaults the other values
// (after get<> was tested)
TYPED_TEST(cartesian_composition_test, custom_assignment)
{
    TypeParam t_d{}; // default to compare

    // first type, default
    TypeParam t1{};

    EXPECT_EQ(get<0>(t1), get<0>(t_d));
    EXPECT_EQ(get<1>(t1), get<1>(t_d));
    EXPECT_NE(get<0>(t1), TestFixture::value_1());
    EXPECT_NE(get<1>(t1), TestFixture::value_2());

    t1 = TestFixture::value_1();

    EXPECT_NE(get<0>(t1), get<0>(t_d));
    EXPECT_EQ(get<1>(t1), get<1>(t_d));
    EXPECT_EQ(get<0>(t1), TestFixture::value_1());
    EXPECT_NE(get<1>(t1), TestFixture::value_2());

    // first type, non-default
    TypeParam t2 = {std::get<4>(TestFixture::values_to_cmp()),
                    std::get<5>(TestFixture::values_to_cmp())};

    EXPECT_NE(get<0>(t2), get<0>(t_d));
    EXPECT_NE(get<1>(t2), get<1>(t_d));
    EXPECT_NE(get<0>(t2), TestFixture::value_1());
    EXPECT_NE(get<1>(t2), TestFixture::value_2());
    EXPECT_EQ(get<0>(t2), std::get<4>(TestFixture::values_to_cmp()));
    EXPECT_EQ(get<1>(t2), std::get<5>(TestFixture::values_to_cmp()));

    t2 = TestFixture::value_1();

    EXPECT_NE(get<0>(t2), get<0>(t_d));
    EXPECT_NE(get<1>(t2), get<1>(t_d));
    EXPECT_EQ(get<0>(t2), TestFixture::value_1());
    EXPECT_NE(get<1>(t2), TestFixture::value_2());
    EXPECT_NE(get<0>(t2), std::get<4>(TestFixture::values_to_cmp()));
    EXPECT_EQ(get<1>(t2), std::get<5>(TestFixture::values_to_cmp()));

    // second type, default
    TypeParam t3{};

    EXPECT_EQ(get<0>(t3), get<0>(t_d));
    EXPECT_EQ(get<1>(t3), get<1>(t_d));
    EXPECT_NE(get<0>(t3), TestFixture::value_1());
    EXPECT_NE(get<1>(t3), TestFixture::value_2());

    t3 = TestFixture::value_2();

    EXPECT_EQ(get<0>(t3), get<0>(t_d));
    EXPECT_NE(get<1>(t3), get<1>(t_d));
    EXPECT_NE(get<0>(t3), TestFixture::value_1());
    EXPECT_EQ(get<1>(t3), TestFixture::value_2());

    // second type, non-default
    TypeParam t4 = {std::get<4>(TestFixture::values_to_cmp()),
                    std::get<5>(TestFixture::values_to_cmp())};

    EXPECT_NE(get<0>(t4), get<0>(t_d));
    EXPECT_NE(get<1>(t4), get<1>(t_d));
    EXPECT_NE(get<0>(t4), TestFixture::value_1());
    EXPECT_NE(get<1>(t4), TestFixture::value_2());
    EXPECT_EQ(get<0>(t4), std::get<4>(TestFixture::values_to_cmp()));
    EXPECT_EQ(get<1>(t4), std::get<5>(TestFixture::values_to_cmp()));

    t4 = TestFixture::value_2();

    EXPECT_NE(get<0>(t4), get<0>(t_d));
    EXPECT_NE(get<1>(t4), get<1>(t_d));
    EXPECT_NE(get<0>(t4), TestFixture::value_1());
    EXPECT_EQ(get<1>(t4), TestFixture::value_2());
    EXPECT_EQ(get<0>(t4), std::get<4>(TestFixture::values_to_cmp()));
    EXPECT_NE(get<1>(t4), std::get<5>(TestFixture::values_to_cmp()));
}

// Custom assignment operator that assigns one type from a subtype and defaults
// the other values
TYPED_TEST(cartesian_composition_test, custom_assignment_subtype)
{
    TypeParam t_d{}; // default to compare

    // first type, default
    TypeParam t1{};

    EXPECT_EQ(get<0>(t1), get<0>(t_d));
    EXPECT_EQ(get<1>(t1), get<1>(t_d));
    EXPECT_NE(get<0>(t1), TestFixture::value_1());
    EXPECT_NE(get<1>(t1), TestFixture::value_2());

    t1 = TestFixture::assignable_to_value_1();

    EXPECT_NE(get<0>(t1), get<0>(t_d));
    EXPECT_EQ(get<1>(t1), get<1>(t_d));
    EXPECT_EQ(get<0>(t1), TestFixture::value_1());
    EXPECT_NE(get<1>(t1), TestFixture::value_2());

    // first type, non-default
    TypeParam t2 = {std::get<4>(TestFixture::values_to_cmp()),
                    std::get<5>(TestFixture::values_to_cmp())};

    EXPECT_NE(get<0>(t2), get<0>(t_d));
    EXPECT_NE(get<1>(t2), get<1>(t_d));
    EXPECT_NE(get<0>(t2), TestFixture::value_1());
    EXPECT_NE(get<1>(t2), TestFixture::value_2());
    EXPECT_EQ(get<0>(t2), std::get<4>(TestFixture::values_to_cmp()));
    EXPECT_EQ(get<1>(t2), std::get<5>(TestFixture::values_to_cmp()));

    t2 = TestFixture::assignable_to_value_1();

    EXPECT_NE(get<0>(t2), get<0>(t_d));
    EXPECT_NE(get<1>(t2), get<1>(t_d));
    EXPECT_EQ(get<0>(t2), TestFixture::value_1());
    EXPECT_NE(get<1>(t2), TestFixture::value_2());
    EXPECT_NE(get<0>(t2), std::get<4>(TestFixture::values_to_cmp()));
    EXPECT_EQ(get<1>(t2), std::get<5>(TestFixture::values_to_cmp()));

    // second type, default
    TypeParam t3{};

    EXPECT_EQ(get<0>(t3), get<0>(t_d));
    EXPECT_EQ(get<1>(t3), get<1>(t_d));
    EXPECT_NE(get<0>(t3), TestFixture::value_1());
    EXPECT_NE(get<1>(t3), TestFixture::value_2());

    t3 = TestFixture::assignable_to_value_2();

    EXPECT_EQ(get<0>(t3), get<0>(t_d));
    EXPECT_NE(get<1>(t3), get<1>(t_d));
    EXPECT_NE(get<0>(t3), TestFixture::value_1());
    EXPECT_EQ(get<1>(t3), TestFixture::value_2());

    // second type, non-default
    TypeParam t4 = {std::get<4>(TestFixture::values_to_cmp()),
                    std::get<5>(TestFixture::values_to_cmp())};

    EXPECT_NE(get<0>(t4), get<0>(t_d));
    EXPECT_NE(get<1>(t4), get<1>(t_d));
    EXPECT_NE(get<0>(t4), TestFixture::value_1());
    EXPECT_NE(get<1>(t4), TestFixture::value_2());
    EXPECT_EQ(get<0>(t4), std::get<4>(TestFixture::values_to_cmp()));
    EXPECT_EQ(get<1>(t4), std::get<5>(TestFixture::values_to_cmp()));

    t4 = TestFixture::assignable_to_value_2();

    EXPECT_NE(get<0>(t4), get<0>(t_d));
    EXPECT_NE(get<1>(t4), get<1>(t_d));
    EXPECT_NE(get<0>(t4), TestFixture::value_1());
    EXPECT_EQ(get<1>(t4), TestFixture::value_2());
    EXPECT_EQ(get<0>(t4), std::get<4>(TestFixture::values_to_cmp()));
    EXPECT_NE(get<1>(t4), std::get<5>(TestFixture::values_to_cmp()));
}

// std::tuple_element
TYPED_TEST(cartesian_composition_test, tuple_element)
{
    using pt = TypeParam;

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, decltype(TestFixture::value_1())>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, decltype(TestFixture::value_2())>);
}

// type deduction
TYPED_TEST(cartesian_composition_test, type_deduce)
{
    TypeParam t0 = TestFixture::instance;
    using pt = decltype(t0);

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, decltype(TestFixture::value_1())>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, decltype(TestFixture::value_2())>);
}

// explicit cast to element
TYPED_TEST(cartesian_composition_test, cast_to_element)
{
    TypeParam t0 = TestFixture::instance;

    auto d = static_cast<decltype(TestFixture::value_1())>(t0);
    auto q = static_cast<decltype(TestFixture::value_2())>(t0);
    static_assert(std::is_same_v<decltype(d), decltype(TestFixture::value_1())>);
    static_assert(std::is_same_v<decltype(q), decltype(TestFixture::value_2())>);

    EXPECT_EQ(d, TestFixture::value_1());
    EXPECT_EQ(q, TestFixture::value_2());
}

// comparison operators
TYPED_TEST(cartesian_composition_test, cmp)
{
    TypeParam t0 = {std::get<2>(TestFixture::values_to_cmp()),
                    std::get<3>(TestFixture::values_to_cmp())};
    TypeParam t1 = {std::get<2>(TestFixture::values_to_cmp()),
                    std::get<1>(TestFixture::values_to_cmp())};
    TypeParam t2 = {std::get<4>(TestFixture::values_to_cmp()),
                    std::get<3>(TestFixture::values_to_cmp())};
    TypeParam t3 = {std::get<0>(TestFixture::values_to_cmp()),
                    std::get<5>(TestFixture::values_to_cmp())};

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

TYPED_TEST(cartesian_composition_test, cmp_to_composite)
{
    // first type
    TypeParam t1 = {std::get<2>(TestFixture::values_to_cmp()),
                    std::get<3>(TestFixture::values_to_cmp())};

    auto [ lt_v1, lt_v2,
           eq_v1, eq_v2,
           gt_v1, gt_v2] = TestFixture::values_to_cmp();

   EXPECT_EQ(t1, eq_v1);
   EXPECT_LE(t1, eq_v1);
   EXPECT_GE(t1, eq_v1);
   EXPECT_LE(t1, gt_v1);
   EXPECT_LT(t1, gt_v1);
   EXPECT_GE(t1, lt_v1);
   EXPECT_GT(t1, lt_v1);

   EXPECT_EQ(eq_v1, t1);
   EXPECT_GE(eq_v1, t1);
   EXPECT_LE(eq_v1, t1);
   EXPECT_GE(gt_v1, t1);
   EXPECT_GT(gt_v1, t1);
   EXPECT_LE(lt_v1, t1);
   EXPECT_LT(lt_v1, t1);

   // second type
    TypeParam t2 = {std::get<2>(TestFixture::values_to_cmp()),
                    std::get<3>(TestFixture::values_to_cmp())};

   EXPECT_EQ(t2, eq_v2);
   EXPECT_LE(t2, eq_v2);
   EXPECT_GE(t2, eq_v2);
   EXPECT_LE(t2, gt_v2);
   EXPECT_LT(t2, gt_v2);
   EXPECT_GE(t2, lt_v2);
   EXPECT_GT(t2, lt_v2);

   EXPECT_EQ(eq_v2, t2);
   EXPECT_GE(eq_v2, t2);
   EXPECT_LE(eq_v2, t2);
   EXPECT_GE(gt_v2, t2);
   EXPECT_GT(gt_v2, t2);
   EXPECT_LE(lt_v2, t2);
   EXPECT_LT(lt_v2, t2);
}

TYPED_TEST(cartesian_composition_test, cmp_to_composite_subtype)
{
    // first type
    TypeParam t0 = {std::get<4>(TestFixture::values_to_cmp()),
                    std::get<5>(TestFixture::values_to_cmp())};
    TypeParam t1 = TestFixture::instance;
    TypeParam t2{};

    EXPECT_EQ(t1, TestFixture::assignable_to_value_1());
    EXPECT_NE(t2, TestFixture::assignable_to_value_1());
    EXPECT_GE(t1, TestFixture::assignable_to_value_1());
    EXPECT_LE(t1, TestFixture::assignable_to_value_1());
    EXPECT_LT(t2, TestFixture::assignable_to_value_1());
    EXPECT_GT(t0, TestFixture::assignable_to_value_1());

    EXPECT_EQ(TestFixture::assignable_to_value_1(), t1);
    EXPECT_NE(TestFixture::assignable_to_value_1(), t0);
    EXPECT_GE(TestFixture::assignable_to_value_1(), t1);
    EXPECT_LE(TestFixture::assignable_to_value_1(), t1);
    EXPECT_LT(TestFixture::assignable_to_value_1(), t0);
    EXPECT_GT(TestFixture::assignable_to_value_1(), t2);

    // second type
    EXPECT_EQ(t1, TestFixture::assignable_to_value_2());
    EXPECT_NE(t2, TestFixture::assignable_to_value_2());
    EXPECT_GE(t1, TestFixture::assignable_to_value_2());
    EXPECT_LE(t1, TestFixture::assignable_to_value_2());
    EXPECT_LT(t2, TestFixture::assignable_to_value_2());
    EXPECT_GT(t0, TestFixture::assignable_to_value_2());

    EXPECT_EQ(TestFixture::assignable_to_value_2(), t1);
    EXPECT_NE(TestFixture::assignable_to_value_2(), t0);
    EXPECT_GE(TestFixture::assignable_to_value_2(), t1);
    EXPECT_LE(TestFixture::assignable_to_value_2(), t1);
    EXPECT_LT(TestFixture::assignable_to_value_2(), t0);
    EXPECT_GT(TestFixture::assignable_to_value_2(), t2);
}
