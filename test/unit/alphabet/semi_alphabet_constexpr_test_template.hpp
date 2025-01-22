// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/concept.hpp>

template <typename t>
using semi_alphabet_constexpr = ::testing::Test;

TYPED_TEST_SUITE_P(semi_alphabet_constexpr);

TYPED_TEST_P(semi_alphabet_constexpr, concept_check)
{
    EXPECT_TRUE(seqan3::detail::constexpr_semialphabet<TypeParam>);
    EXPECT_TRUE(seqan3::detail::constexpr_semialphabet<TypeParam &>);

    EXPECT_TRUE(seqan3::detail::constexpr_semialphabet<TypeParam const>);
    EXPECT_TRUE(seqan3::detail::constexpr_semialphabet<TypeParam const &>);

    EXPECT_TRUE(seqan3::detail::writable_constexpr_semialphabet<TypeParam>);
    EXPECT_TRUE(seqan3::detail::writable_constexpr_semialphabet<TypeParam &>);

    EXPECT_FALSE(seqan3::detail::writable_constexpr_semialphabet<TypeParam const>);
    EXPECT_FALSE(seqan3::detail::writable_constexpr_semialphabet<TypeParam const &>);
}

TYPED_TEST_P(semi_alphabet_constexpr, default_value_constructor)
{
    [[maybe_unused]] constexpr TypeParam t0{};
}

TYPED_TEST_P(semi_alphabet_constexpr, assign_rank)
{
    constexpr size_t rank = 1 % seqan3::alphabet_size<TypeParam>;
    [[maybe_unused]] constexpr TypeParam t0{seqan3::assign_rank_to(rank, TypeParam{})};
}

TYPED_TEST_P(semi_alphabet_constexpr, to_rank)
{
    constexpr size_t rank = 1 % seqan3::alphabet_size<TypeParam>;
    constexpr TypeParam t0{seqan3::assign_rank_to(rank, TypeParam{})};
    constexpr bool b = (seqan3::to_rank(t0) == rank);
    EXPECT_TRUE(b);
}

TYPED_TEST_P(semi_alphabet_constexpr, copy_constructor)
{
    constexpr seqan3::alphabet_rank_t<TypeParam> rank = 1 % seqan3::alphabet_size<TypeParam>;
    constexpr TypeParam t1{seqan3::assign_rank_to(rank, TypeParam{})};

    constexpr TypeParam t2{t1};
    constexpr TypeParam t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

TYPED_TEST_P(semi_alphabet_constexpr, move_constructor)
{
    constexpr seqan3::alphabet_rank_t<TypeParam> rank = 1 % seqan3::alphabet_size<TypeParam>;
    constexpr TypeParam t0{seqan3::assign_rank_to(rank, TypeParam{})};
    constexpr TypeParam t1{t0};

    constexpr TypeParam t2{std::move(t1)};
    constexpr TypeParam t3(std::move(t2));
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t3, t0);
}

TYPED_TEST_P(semi_alphabet_constexpr, copy_assignment)
{
    constexpr size_t rank = 1 % seqan3::alphabet_size<TypeParam>;
    constexpr TypeParam t0{seqan3::assign_rank_to(rank, TypeParam{})};
    // constexpr context:
    constexpr TypeParam t3 = [&]() constexpr
    {
        TypeParam t1{seqan3::assign_rank_to(rank, TypeParam{})};
        TypeParam t2{};
        t2 = t1;

        return t2;
    }();
    EXPECT_EQ(t3, t0);
}

TYPED_TEST_P(semi_alphabet_constexpr, move_assignment)
{
    constexpr size_t rank = 1 % seqan3::alphabet_size<TypeParam>;
    constexpr TypeParam t0{seqan3::assign_rank_to(rank, TypeParam{})};
    // constexpr context:
    constexpr TypeParam t3 = [&]() constexpr
    {
        TypeParam t1{seqan3::assign_rank_to(rank, TypeParam{})};
        TypeParam t2{};
        t2 = std::move(t1);

        return t2;
    }();
    EXPECT_EQ(t3, t0);
}

TYPED_TEST_P(semi_alphabet_constexpr, comparison_operators)
{
    if constexpr (seqan3::alphabet_size<TypeParam> == 1)
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
        constexpr TypeParam t0{seqan3::assign_rank_to(0, TypeParam{})};
        constexpr TypeParam t1{seqan3::assign_rank_to(1, TypeParam{})};
        constexpr bool b1 = (t0 < t1);
        constexpr bool b2 = (t0 <= t1);
        constexpr bool b3 = (t1 <= t1);
        constexpr bool b4 = (t1 == t1);
        constexpr bool b5 = (t1 >= t1);
        constexpr bool b6 = (t1 >= t0);
        constexpr bool b7 = (t1 > t0);
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

REGISTER_TYPED_TEST_SUITE_P(semi_alphabet_constexpr,
                            concept_check,
                            default_value_constructor,
                            assign_rank,
                            to_rank,
                            copy_constructor,
                            move_constructor,
                            copy_assignment,
                            move_assignment,
                            comparison_operators);
