// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/exception.hpp>
#include <seqan3/test/pretty_printing.hpp>
#include <seqan3/utility/concept.hpp>

template <typename T>
using semi_alphabet_test = ::testing::Test;

constexpr size_t maximum_iterations = 65536u;

TYPED_TEST_SUITE_P(semi_alphabet_test);

TYPED_TEST_P(semi_alphabet_test, concept_check)
{
    EXPECT_TRUE(seqan3::semialphabet<TypeParam>);
    EXPECT_TRUE(seqan3::semialphabet<TypeParam &>);
    EXPECT_TRUE(seqan3::semialphabet<TypeParam const>);
    EXPECT_TRUE(seqan3::semialphabet<TypeParam const &>);

    EXPECT_TRUE(seqan3::writable_semialphabet<TypeParam>);
    EXPECT_TRUE(seqan3::writable_semialphabet<TypeParam &>);
    EXPECT_FALSE(seqan3::writable_semialphabet<TypeParam const>);
    EXPECT_FALSE(seqan3::writable_semialphabet<TypeParam const &>);
}

TYPED_TEST_P(semi_alphabet_test, type_properties)
{
    // It is highly recommended that non-reference types that model this concept, also model:
    EXPECT_TRUE((std::regular<TypeParam>));
    EXPECT_TRUE((seqan3::trivially_copyable<TypeParam>));
    EXPECT_TRUE((seqan3::standard_layout<TypeParam>));
}

TYPED_TEST_P(semi_alphabet_test, alphabet_size)
{
    EXPECT_GT(seqan3::alphabet_size<TypeParam>, 0u);
}

TYPED_TEST_P(semi_alphabet_test, default_value_constructor)
{
    EXPECT_TRUE((std::is_nothrow_default_constructible_v<TypeParam>));
}

TYPED_TEST_P(semi_alphabet_test, assign_rank_to)
{
    // this double checks the value initialisation
    EXPECT_EQ((seqan3::assign_rank_to(0, TypeParam{})), TypeParam{});

    TypeParam t0;
    for (size_t i = 0u; i < seqan3::alphabet_size<TypeParam> && i < maximum_iterations; ++i)
        seqan3::assign_rank_to(i, t0);

    EXPECT_TRUE((std::is_same_v<decltype(seqan3::assign_rank_to(0, t0)), TypeParam &>));
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::assign_rank_to(0, TypeParam{})), TypeParam>));
}

TYPED_TEST_P(semi_alphabet_test, to_rank)
{
    // this double checks the value initialisation
    EXPECT_EQ(seqan3::to_rank(TypeParam{}), 0u);

    TypeParam t0;
    for (size_t i = 0; i < seqan3::alphabet_size<TypeParam> && i < maximum_iterations; ++i)
        EXPECT_EQ((seqan3::to_rank(seqan3::assign_rank_to(i, t0))), i);

    EXPECT_TRUE((std::is_same_v<decltype(seqan3::to_rank(t0)), seqan3::alphabet_rank_t<TypeParam>>));
}

TYPED_TEST_P(semi_alphabet_test, copy_constructor)
{
    // the module operation ensures that the result is within the valid rank range;
    // it will be in the most cases 1 except for alphabets like seqan3::gap where it will be 0
    constexpr seqan3::alphabet_rank_t<TypeParam> rank = 1 % seqan3::alphabet_size<TypeParam>;
    TypeParam t1;
    seqan3::assign_rank_to(rank, t1);
    TypeParam t2{t1};
    TypeParam t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

TYPED_TEST_P(semi_alphabet_test, move_constructor)
{
    constexpr seqan3::alphabet_rank_t<TypeParam> rank = 1 % seqan3::alphabet_size<TypeParam>;
    TypeParam t0;
    seqan3::assign_rank_to(rank, t0);
    TypeParam t1{t0};

    TypeParam t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    TypeParam t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

TYPED_TEST_P(semi_alphabet_test, copy_assignment)
{
    constexpr seqan3::alphabet_rank_t<TypeParam> rank = 1 % seqan3::alphabet_size<TypeParam>;
    TypeParam t1;
    seqan3::assign_rank_to(rank, t1);
    TypeParam t2;
    t2 = t1;
    EXPECT_EQ(t1, t2);
}

TYPED_TEST_P(semi_alphabet_test, move_assignment)
{
    constexpr seqan3::alphabet_rank_t<TypeParam> rank = 1 % seqan3::alphabet_size<TypeParam>;
    TypeParam t0;
    seqan3::assign_rank_to(rank, t0);
    TypeParam t1{t0};
    TypeParam t2;
    TypeParam t3;
    t2 = std::move(t1);
    EXPECT_EQ(t2, t0);
    t3 = std::move(t2);
    EXPECT_EQ(t3, t0);
}

TYPED_TEST_P(semi_alphabet_test, swap)
{
    constexpr seqan3::alphabet_rank_t<TypeParam> rank = 1 % seqan3::alphabet_size<TypeParam>;
    TypeParam t0;
    seqan3::assign_rank_to(rank, t0);
    TypeParam t1{t0};
    TypeParam t2{};
    TypeParam t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

TYPED_TEST_P(semi_alphabet_test, comparison_operators)
{
    TypeParam t0{};
    TypeParam t1{};

    seqan3::assign_rank_to(0, t0);
    seqan3::assign_rank_to(1 % seqan3::alphabet_size<TypeParam>, t1);

    EXPECT_EQ(t0, t0);
    EXPECT_LE(t0, t1);
    EXPECT_LE(t1, t1);
    EXPECT_EQ(t1, t1);
    EXPECT_GE(t1, t1);
    EXPECT_GE(t1, t0);

    if constexpr (seqan3::alphabet_size<TypeParam> == 1)
    {
        EXPECT_EQ(t0, t1);
    }
    else
    {
        EXPECT_LT(t0, t1);
        EXPECT_NE(t0, t1);
        EXPECT_GT(t1, t0);
    }
}

REGISTER_TYPED_TEST_SUITE_P(semi_alphabet_test,
                            concept_check,
                            type_properties,
                            alphabet_size,
                            default_value_constructor,
                            assign_rank_to,
                            to_rank,
                            copy_constructor,
                            move_constructor,
                            copy_assignment,
                            move_assignment,
                            swap,
                            comparison_operators);
