// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/debug_stream.hpp>

#if SEQAN3_WITH_CEREAL
#include <seqan3/test/tmp_filename.hpp>

#include <fstream>

#include <cereal/archives/binary.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/xml.hpp>
#include <cereal/types/vector.hpp>
#endif // SEQAN3_WITH_CEREAL

using namespace seqan3;

template <typename T>
using alphabet_constexpr = ::testing::Test;

TYPED_TEST_CASE_P(alphabet_constexpr);

TYPED_TEST_P(alphabet_constexpr, concept_check)
{
    EXPECT_TRUE(detail::ConstexprAlphabet<TypeParam   >);
    EXPECT_TRUE(detail::ConstexprAlphabet<TypeParam & >);

    EXPECT_TRUE(detail::ConstexprAlphabet<TypeParam const   >);
    EXPECT_TRUE(detail::ConstexprAlphabet<TypeParam const & >);

    EXPECT_TRUE(detail::WritableConstexprAlphabet<TypeParam   >);
    EXPECT_TRUE(detail::WritableConstexprAlphabet<TypeParam & >);

    EXPECT_FALSE(detail::WritableConstexprAlphabet<TypeParam const   >);
    EXPECT_FALSE(detail::WritableConstexprAlphabet<TypeParam const & >);
}

TYPED_TEST_P(alphabet_constexpr, default_value_constructor)
{
    [[maybe_unused]] constexpr TypeParam t0{};
}

TYPED_TEST_P(alphabet_constexpr, copy_constructor)
{
    constexpr TypeParam t1{};
    constexpr TypeParam t2{t1};
    constexpr TypeParam t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

TYPED_TEST_P(alphabet_constexpr, move_constructor)
{
    constexpr TypeParam t0{};
    constexpr TypeParam t1{t0};

    constexpr TypeParam t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    constexpr TypeParam t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

TYPED_TEST_P(alphabet_constexpr, global_assign_rank)
{
    constexpr size_t rank = 1 % alphabet_size<TypeParam>;
    [[maybe_unused]] constexpr TypeParam t0{assign_rank_to(rank, TypeParam{})};
}

TYPED_TEST_P(alphabet_constexpr, global_to_rank)
{
    constexpr size_t rank = 1 % alphabet_size<TypeParam>;
    constexpr TypeParam t0{assign_rank_to(rank, TypeParam{})};
    constexpr bool b = (to_rank(t0) == rank);
    EXPECT_TRUE(b);
}

TYPED_TEST_P(alphabet_constexpr, copy_assignment)
{
    constexpr size_t rank = 1 % alphabet_size<TypeParam>;
    constexpr TypeParam t0{assign_rank_to(rank, TypeParam{})};
    // constexpr context:
    constexpr TypeParam t3 = [&] () constexpr
    {
        TypeParam t1{assign_rank_to(rank, TypeParam{})};
        TypeParam t2{};
        t2 = t1;

        return t2;
    }();
    EXPECT_EQ(t3, t0);
}

TYPED_TEST_P(alphabet_constexpr, move_assignment)
{
    constexpr size_t rank = 1 % alphabet_size<TypeParam>;
    constexpr TypeParam t0{assign_rank_to(rank, TypeParam{})};
    // constexpr context:
    constexpr TypeParam t3 = [&] () constexpr
    {
        TypeParam t1{assign_rank_to(rank, TypeParam{})};
        TypeParam t2{};
        t2 = std::move(t1);

        return t2;
    }();
    EXPECT_EQ(t3, t0);
}

TYPED_TEST_P(alphabet_constexpr, global_assign_char)
{
    [[maybe_unused]] constexpr TypeParam t0{assign_char_to('A', TypeParam{})};
}

TYPED_TEST_P(alphabet_constexpr, global_to_char)
{
    constexpr TypeParam t0{};
    [[maybe_unused]] constexpr alphabet_char_t<TypeParam> c = to_char(t0);
}

TYPED_TEST_P(alphabet_constexpr, comparison_operators)
{
    if constexpr (alphabet_size<TypeParam> == 1)
    {
        constexpr TypeParam t0{};
        constexpr TypeParam t1{};
        constexpr bool b2 = (t0 <= t1);
        constexpr bool b3 = (t1 <= t1);
        constexpr bool b4 = (t1 == t1);
        constexpr bool b5 = (t1 >= t1);
        constexpr bool b6 = (t1 >= t0);

        EXPECT_TRUE(b2);
        EXPECT_TRUE(b3);
        EXPECT_TRUE(b4);
        EXPECT_TRUE(b5);
        EXPECT_TRUE(b6);
    }
    else
    {
        constexpr TypeParam t0{assign_rank_to(0, TypeParam{})};
        constexpr TypeParam t1{assign_rank_to(1, TypeParam{})};
        constexpr bool b1 = (t0 <  t1);
        constexpr bool b2 = (t0 <= t1);
        constexpr bool b3 = (t1 <= t1);
        constexpr bool b4 = (t1 == t1);
        constexpr bool b5 = (t1 >= t1);
        constexpr bool b6 = (t1 >= t0);
        constexpr bool b7 = (t1 >  t0);
        constexpr bool b8 = (t0 != t1);

        EXPECT_TRUE(b1);
        EXPECT_TRUE(b2);
        EXPECT_TRUE(b3);
        EXPECT_TRUE(b4);
        EXPECT_TRUE(b5);
        EXPECT_TRUE(b6);
        EXPECT_TRUE(b7);
        EXPECT_TRUE(b8);
    }
}

REGISTER_TYPED_TEST_CASE_P(alphabet_constexpr, concept_check, default_value_constructor, copy_constructor,
    move_constructor, global_assign_rank, global_to_rank, copy_assignment, move_assignment, global_assign_char,
    global_to_char, comparison_operators);
