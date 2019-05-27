// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/concept/tuple.hpp>

using namespace seqan3;

template <typename type1, typename type2>
struct test_composite : public alphabet_tuple_base<test_composite<type1, type2>, type1, type2>
{
    using base_t = alphabet_tuple_base<test_composite<type1, type2>, type1, type2>;
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
class alphabet_tuple_base_test : public ::testing::Test {};

template <>
class alphabet_tuple_base_test<test_composite<dna4, dna5>> : public ::testing::Test
{
public:
    using T = test_composite<dna4, dna5>;

    T instance = T{value_1(), value_2()};
    T zero_instance = T{decltype(value_1()){}, decltype(value_2()){}};
    size_t tup_size{2};

    // test_composite<dna4, dna5>
    // -------------------------------------------------------------------------
    dna4 value_1()
    {
        return 'G'_dna4;
    }
    rna4 assignable_to_value_1()
    {
        return 'G'_rna4;
    }
    dna5 value_2()
    {
        return 'G'_dna5;
    }
    rna5 assignable_to_value_2()
    {
        return 'G'_rna5;
    }
    auto values_to_cmp()
    {
        return std::make_tuple(/*low */'A'_dna4, 'A'_dna5,
                               /*mid */'C'_dna4, 'C'_dna5,
                               /*high*/'T'_dna4, 'T'_dna5);
    }
};

template <>
class alphabet_tuple_base_test<qualified<dna4, phred42>>: public ::testing::Test
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
        return 'G'_dna4;
    }
    rna4 assignable_to_value_1()
    {
        return 'G'_rna4;
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
        return std::make_tuple(/*low */'A'_dna4, phred42{1},
                               /*mid */'C'_dna4, phred42{4},
                               /*high*/'T'_dna4, phred42{9});
    }
};

template <>
class alphabet_tuple_base_test<structured_rna<rna4, dot_bracket3>>: public ::testing::Test
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
        return 'G'_rna4;
    }
    dna4 assignable_to_value_1()
    {
        return 'G'_dna4;
    }
    dot_bracket3 value_2()
    {
        return '('_db3;
    }
    dot_bracket3 assignable_to_value_2()
    {
        return '('_db3; // replace if assignable subtype becomes available
    }
    auto values_to_cmp()
    {
        return std::make_tuple(/*low */'A'_rna4, '.'_db3,
                               /*mid */'C'_rna4, '('_db3,
                               /*high*/'T'_rna4, ')'_db3);
    }
};

template <>
class alphabet_tuple_base_test<structured_aa<aa27, dssp9>>: public ::testing::Test
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
        return 'K'_aa27;
    }
    aa27 assignable_to_value_1()
    {
        return 'K'_aa27; // replace if assignable subtype becomes available
    }
    dssp9 value_2()
    {
        return 'I'_dssp9;
    }
    dssp9 assignable_to_value_2()
    {
        return 'I'_dssp9; // replace if assignable subtype becomes available
    }
    auto values_to_cmp()
    {
        return std::make_tuple(/*low */'A'_aa27, 'H'_dssp9,
                               /*mid */'P'_aa27, 'I'_dssp9,
                               /*high*/'Z'_aa27, 'X'_dssp9);
    }
};

using composite_types = ::testing::Types<test_composite<dna4, dna5>,
                                           structured_rna<rna4, dot_bracket3>,
                                           structured_aa<aa27, dssp9>,
                                           qualified<dna4, phred42>>;

TYPED_TEST_CASE(alphabet_tuple_base_test, composite_types);

TYPED_TEST(alphabet_tuple_base_test, concept_check)
{
    EXPECT_TRUE(TupleLike<TypeParam>);
}

// default/zero construction
TYPED_TEST(alphabet_tuple_base_test, ctr)
{
    [[maybe_unused]] TypeParam t1{};
    EXPECT_EQ(std::tuple_size<TypeParam>::value, TestFixture::tup_size);
}

// initialiser-list initialization
TYPED_TEST(alphabet_tuple_base_test, aggr)
{
    TypeParam t1{};
    TypeParam t2 = TestFixture::instance; // test in fixture to be type independent

    EXPECT_NE(t1, t2);
}

// copy assignment
TYPED_TEST(alphabet_tuple_base_test, cp_assgn)
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
TYPED_TEST(alphabet_tuple_base_test, zro)
{
    TypeParam t1 = TestFixture::zero_instance;
    TypeParam t2{};

    EXPECT_EQ(t1, t2);
}

// copy construction
TYPED_TEST(alphabet_tuple_base_test, cp_ctr)
{
    TypeParam t1 = TestFixture::instance;
    TypeParam t2{t1};
    TypeParam t3(t1);

    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move construction
TYPED_TEST(alphabet_tuple_base_test, mv_ctr)
{
    TypeParam t0 = TestFixture::instance;
    TypeParam t1 = TestFixture::instance;
    TypeParam t2{std::move(t1)};

    EXPECT_EQ(t2, t0);

    TypeParam t3(std::move(t2));

    EXPECT_EQ(t3, t0);
}

// move assignment
TYPED_TEST(alphabet_tuple_base_test, mv_assgn)
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
TYPED_TEST(alphabet_tuple_base_test, swap)
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
TYPED_TEST(alphabet_tuple_base_test, get_i)
{
    TypeParam t0 = TestFixture::instance;

//     static_assert(std::is_same_v<decltype(get<0>(t0)), decltype(TestFixture::value_1()) &>);
//     static_assert(std::is_same_v<decltype(get<1>(t0)), decltype(TestFixture::value_2()) &>);

    EXPECT_EQ(get<0>(t0), TestFixture::value_1());
    EXPECT_EQ(get<1>(t0), TestFixture::value_2());
}

// structured bindings
TYPED_TEST(alphabet_tuple_base_test, struct_binding)
{
  TypeParam t0 = TestFixture::instance;
  auto [ i, l ] = t0;

  static_assert(std::is_same_v<decltype(i), decltype(TestFixture::value_1())>);
  static_assert(std::is_same_v<decltype(l), decltype(TestFixture::value_2())>);

  EXPECT_EQ(i, TestFixture::value_1());
  EXPECT_EQ(l, TestFixture::value_2());
}

// get<type>
TYPED_TEST(alphabet_tuple_base_test, get_type)
{
    TypeParam t0 = TestFixture::instance;

    EXPECT_EQ(get<decltype(TestFixture::value_1())>(t0), TestFixture::value_1());
    EXPECT_EQ(get<decltype(TestFixture::value_2())>(t0), TestFixture::value_2());
}

// Custom constructor that assigns one type and defaults the other values
// (after get<> was tested)
TYPED_TEST(alphabet_tuple_base_test, custom_ctr)
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
TYPED_TEST(alphabet_tuple_base_test, custom_ctr_subtype)
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
TYPED_TEST(alphabet_tuple_base_test, custom_assignment)
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
TYPED_TEST(alphabet_tuple_base_test, custom_assignment_subtype)
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
TYPED_TEST(alphabet_tuple_base_test, tuple_element)
{
    using pt = TypeParam;

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, decltype(TestFixture::value_1())>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, decltype(TestFixture::value_2())>);
}

// type deduction
TYPED_TEST(alphabet_tuple_base_test, type_deduce)
{
    TypeParam t0 = TestFixture::instance;
    using pt = decltype(t0);

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, decltype(TestFixture::value_1())>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, decltype(TestFixture::value_2())>);
}

// explicit cast to element
TYPED_TEST(alphabet_tuple_base_test, cast_to_element)
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
TYPED_TEST(alphabet_tuple_base_test, cmp)
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

TYPED_TEST(alphabet_tuple_base_test, cmp_to_composite)
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

TYPED_TEST(alphabet_tuple_base_test, cmp_to_composite_subtype)
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
