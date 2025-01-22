// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>

#include <seqan3/alphabet/adaptation/char.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

using char_types = ::testing::Types<char, char16_t, char32_t, wchar_t>;

INSTANTIATE_TYPED_TEST_SUITE_P(char_adaptation, alphabet, char_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(char_adaptation, semi_alphabet_test, char_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(char_adaptation, alphabet_constexpr, char_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(char_adaptation, semi_alphabet_constexpr, char_types, );

template <typename T>
using char_adaptation = ::testing::Test;

TYPED_TEST_SUITE(char_adaptation, char_types, );

TYPED_TEST(char_adaptation, type_properties)
{
    EXPECT_TRUE((std::is_trivially_copyable_v<TypeParam>));
    EXPECT_TRUE((std::is_trivially_default_constructible_v<TypeParam>));
    EXPECT_TRUE((seqan3::trivial<TypeParam>));
}

TYPED_TEST(char_adaptation, alphabet_char_t)
{
    EXPECT_TRUE((std::is_same_v<seqan3::alphabet_char_t<TypeParam>, TypeParam>));
    EXPECT_TRUE((std::is_same_v<seqan3::alphabet_char_t<TypeParam &>, TypeParam>));
    EXPECT_TRUE((std::is_same_v<seqan3::alphabet_char_t<TypeParam &&>, TypeParam>));
}

TYPED_TEST(char_adaptation, to_char)
{
    TypeParam l{'A'};
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::to_char(l)), seqan3::alphabet_char_t<TypeParam>>));
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::to_char(TypeParam{'A'})), seqan3::alphabet_char_t<TypeParam>>));
    EXPECT_EQ(seqan3::to_char(TypeParam{'A'}), l);
}

TYPED_TEST(char_adaptation, assign_char_to)
{
    TypeParam l{'A'};
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::assign_char_to('A', l)), seqan3::alphabet_char_t<TypeParam> &>));
    EXPECT_TRUE(
        (std::is_same_v<decltype(seqan3::assign_char_to('A', TypeParam{'A'})), seqan3::alphabet_char_t<TypeParam>>));
    EXPECT_EQ((seqan3::assign_char_to('A', TypeParam{'C'})), l);
    EXPECT_EQ((seqan3::assign_char_to('C', l)), TypeParam{'C'});
}

TYPED_TEST(char_adaptation, assign_char_strictly_to)
{
    TypeParam l{'A'};
    EXPECT_TRUE(
        (std::is_same_v<decltype(seqan3::assign_char_strictly_to('A', l)), seqan3::alphabet_char_t<TypeParam> &>));
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::assign_char_strictly_to('A', TypeParam{'A'})),
                                seqan3::alphabet_char_t<TypeParam>>));
    EXPECT_EQ((seqan3::assign_char_strictly_to('A', TypeParam{'C'})), l);
    EXPECT_EQ((seqan3::assign_char_strictly_to('C', l)), TypeParam{'C'});
}

TYPED_TEST(char_adaptation, alphabet_rank_t)
{
    EXPECT_TRUE((std::is_integral_v<seqan3::alphabet_rank_t<TypeParam>>));
    EXPECT_TRUE((std::is_unsigned_v<seqan3::alphabet_rank_t<TypeParam>>));
    EXPECT_GE(sizeof(seqan3::alphabet_rank_t<TypeParam>), sizeof(TypeParam));
}

TYPED_TEST(char_adaptation, to_rank)
{
    TypeParam l{'A'};
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::to_rank(l)), seqan3::alphabet_rank_t<TypeParam>>));
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::to_rank(TypeParam{'A'})), seqan3::alphabet_rank_t<TypeParam>>));

    unsigned char cmp{'A'};
    EXPECT_EQ(seqan3::to_rank(TypeParam{65}), cmp);
}

TYPED_TEST(char_adaptation, assign_rank_to)
{
    TypeParam l{'A'};
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::assign_rank_to(65, l)), seqan3::alphabet_char_t<TypeParam> &>));
    EXPECT_TRUE(
        (std::is_same_v<decltype(seqan3::assign_rank_to(65, TypeParam{'A'})), seqan3::alphabet_char_t<TypeParam>>));
    EXPECT_EQ((seqan3::assign_rank_to(65, TypeParam{'C'})), l);
    EXPECT_EQ((seqan3::assign_rank_to(67, l)), TypeParam{'C'});
}

TYPED_TEST(char_adaptation, alphabet_size)
{
    EXPECT_EQ(seqan3::alphabet_size<TypeParam>,
              static_cast<size_t>(std::numeric_limits<TypeParam>::max()) + 1
                  - std::numeric_limits<TypeParam>::lowest());
}
