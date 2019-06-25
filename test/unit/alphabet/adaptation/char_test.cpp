// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/adaptation/all.hpp>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"

using namespace seqan3;

// char32_t and wchar_t, too slow
using fast_char_types = ::testing::Types<char, char16_t, char32_t, wchar_t>;

INSTANTIATE_TYPED_TEST_CASE_P(char_adaptation, alphabet, fast_char_types);
INSTANTIATE_TYPED_TEST_CASE_P(char_adaptation, alphabet_constexpr, fast_char_types);

template <typename T>
using char_adaptation = ::testing::Test;

using char_types = ::testing::Types<char, char16_t, char32_t, wchar_t>;

TYPED_TEST_CASE(char_adaptation, char_types);

TYPED_TEST(char_adaptation, type_properties)
{
    EXPECT_TRUE((std::is_trivially_copyable_v<TypeParam>));
    EXPECT_TRUE((std::is_trivially_default_constructible_v<TypeParam>));
    EXPECT_TRUE((std::is_trivial_v<TypeParam>));
}

TYPED_TEST(char_adaptation, alphabet_char_t)
{
    EXPECT_TRUE((std::is_same_v<alphabet_char_t<TypeParam   >, TypeParam>));
    EXPECT_TRUE((std::is_same_v<alphabet_char_t<TypeParam & >, TypeParam>));
    EXPECT_TRUE((std::is_same_v<alphabet_char_t<TypeParam &&>, TypeParam>));
}

TYPED_TEST(char_adaptation, to_char)
{
    TypeParam l{'A'};
    EXPECT_TRUE((std::is_same_v<decltype(to_char(l)),              alphabet_char_t<TypeParam>>));
    EXPECT_TRUE((std::is_same_v<decltype(to_char(TypeParam{'A'})), alphabet_char_t<TypeParam>>));
    EXPECT_EQ(to_char(TypeParam{'A'}), l);
}

TYPED_TEST(char_adaptation, assign_char_to)
{
    TypeParam l{'A'};
    EXPECT_TRUE((std::is_same_v<decltype(assign_char_to('A', l             )), alphabet_char_t<TypeParam> &>));
    EXPECT_TRUE((std::is_same_v<decltype(assign_char_to('A', TypeParam{'A'})), alphabet_char_t<TypeParam>  >));
    EXPECT_EQ((assign_char_to('A', TypeParam{'C'})), l);
    EXPECT_EQ((assign_char_to('C', l             )), TypeParam{'C'});
}

TYPED_TEST(char_adaptation, assign_char_strictly_to)
{
    TypeParam l{'A'};
    EXPECT_TRUE((std::is_same_v<decltype(assign_char_strictly_to('A', l             )), alphabet_char_t<TypeParam> &>));
    EXPECT_TRUE((std::is_same_v<decltype(assign_char_strictly_to('A', TypeParam{'A'})), alphabet_char_t<TypeParam>  >));
    EXPECT_EQ((assign_char_strictly_to('A', TypeParam{'C'})), l);
    EXPECT_EQ((assign_char_strictly_to('C', l             )), TypeParam{'C'});
}

TYPED_TEST(char_adaptation, alphabet_rank_t)
{
    EXPECT_TRUE((std::is_integral_v<alphabet_rank_t<TypeParam>>));
    EXPECT_TRUE((std::is_unsigned_v<alphabet_rank_t<TypeParam>>));
    EXPECT_GE(sizeof(alphabet_rank_t<TypeParam>), sizeof(TypeParam));
}

TYPED_TEST(char_adaptation, to_rank)
{
    TypeParam l{'A'};
    EXPECT_TRUE((std::is_same_v<decltype(to_rank(l)),              alphabet_rank_t<TypeParam>>));
    EXPECT_TRUE((std::is_same_v<decltype(to_rank(TypeParam{'A'})), alphabet_rank_t<TypeParam>>));

    unsigned char cmp{'A'};
    EXPECT_EQ(to_rank(TypeParam{65}), cmp);
}

TYPED_TEST(char_adaptation, assign_rank_to)
{
    TypeParam l{'A'};
    EXPECT_TRUE((std::is_same_v<decltype(assign_rank_to(65, l)), alphabet_char_t<TypeParam> &>));
    EXPECT_TRUE((std::is_same_v<decltype(assign_rank_to(65, TypeParam{'A'})), alphabet_char_t<TypeParam>  >));
    EXPECT_EQ((assign_rank_to(65, TypeParam{'C'})), l);
    EXPECT_EQ((assign_rank_to(67, l)), TypeParam{'C'});
}

TYPED_TEST(char_adaptation, alphabet_size)
{
    EXPECT_EQ(alphabet_size<TypeParam>,
        static_cast<size_t>(std::numeric_limits<TypeParam>::max()) + 1 - std::numeric_limits<TypeParam>::lowest());
}
