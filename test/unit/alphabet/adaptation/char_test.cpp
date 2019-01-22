// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/adaptation/all.hpp>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"

using namespace seqan3;

// char32_t and wchar_t, too slow
using fast_char_types = ::testing::Types<char, char16_t/*, char32_t, wchar_t*/>;

INSTANTIATE_TYPED_TEST_CASE_P(char_adaptation, alphabet, fast_char_types);
INSTANTIATE_TYPED_TEST_CASE_P(char_adaptation, alphabet_constexpr, fast_char_types);

template <typename T>
class char_adaptation : public ::testing::Test
{};

using char_types = ::testing::Types<char, char16_t, char32_t, wchar_t>;

TYPED_TEST_CASE(char_adaptation, char_types);

TYPED_TEST(char_adaptation, concept_check)
{
    EXPECT_TRUE(char_adaptation_concept<TypeParam>);
    // NOTE: Using intermediate concept notation with forwarding references cause the concept type
    // to hold a reference.
    EXPECT_TRUE(char_adaptation_concept<TypeParam &>);
    EXPECT_TRUE(char_adaptation_concept<TypeParam &&>);
}

TYPED_TEST(char_adaptation, underlying_char_t)
{
    EXPECT_TRUE((std::is_same_v<underlying_char_t<TypeParam   >, TypeParam>));
    EXPECT_TRUE((std::is_same_v<underlying_char_t<TypeParam & >, TypeParam>));
    EXPECT_TRUE((std::is_same_v<underlying_char_t<TypeParam &&>, TypeParam>));
}

TYPED_TEST(char_adaptation, to_char)
{
    TypeParam l{'A'};
    EXPECT_TRUE((std::is_same_v<decltype(to_char(l)),              underlying_char_t<TypeParam>>));
    EXPECT_TRUE((std::is_same_v<decltype(to_char(TypeParam{'A'})), underlying_char_t<TypeParam>>));
    EXPECT_EQ(to_char(TypeParam{'A'}), l);
}

TYPED_TEST(char_adaptation, assign_char)
{
    TypeParam l{'A'};
    EXPECT_TRUE((std::is_same_v<decltype(assign_char(l,              'A')), underlying_char_t<TypeParam> &>));
    EXPECT_TRUE((std::is_same_v<decltype(assign_char(TypeParam{'A'}, 'A')), underlying_char_t<TypeParam>  >));
    EXPECT_EQ((assign_char(TypeParam{'C'}, 'A')), l);
    EXPECT_EQ((assign_char(l,              'C')), TypeParam{'C'});
}

TYPED_TEST(char_adaptation, assign_char_strict)
{
    TypeParam l{'A'};
    EXPECT_TRUE((std::is_same_v<decltype(assign_char_strict(l,              'A')), underlying_char_t<TypeParam> &>));
    EXPECT_TRUE((std::is_same_v<decltype(assign_char_strict(TypeParam{'A'}, 'A')), underlying_char_t<TypeParam>  >));
    EXPECT_EQ((assign_char_strict(TypeParam{'C'}, 'A')), l);
    EXPECT_EQ((assign_char_strict(l,              'C')), TypeParam{'C'});
}

TYPED_TEST(char_adaptation, underlying_rank_t)
{
    EXPECT_TRUE((std::is_integral_v<underlying_rank_t<TypeParam>>));
    EXPECT_TRUE((std::is_unsigned_v<underlying_rank_t<TypeParam>>));
    EXPECT_GE(sizeof(underlying_rank_t<TypeParam>), sizeof(TypeParam));
}

TYPED_TEST(char_adaptation, to_rank)
{
    TypeParam l{'A'};
    EXPECT_TRUE((std::is_same_v<decltype(to_rank(l)),              underlying_rank_t<TypeParam>>));
    EXPECT_TRUE((std::is_same_v<decltype(to_rank(TypeParam{'A'})), underlying_rank_t<TypeParam>>));

    unsigned char cmp{'A'};
    EXPECT_EQ(to_rank(TypeParam{65}), cmp);
}

TYPED_TEST(char_adaptation, assign_rank)
{
    TypeParam l{'A'};
    EXPECT_TRUE((std::is_same_v<decltype(assign_rank(l,              65)), underlying_char_t<TypeParam> &>));
    EXPECT_TRUE((std::is_same_v<decltype(assign_rank(TypeParam{'A'}, 65)), underlying_char_t<TypeParam>  >));
    EXPECT_EQ((assign_rank(TypeParam{'C'}, 65)), l);
    EXPECT_EQ((assign_rank(l,              67)), TypeParam{'C'});
}

TYPED_TEST(char_adaptation, alphabet_size_v)
{
    EXPECT_EQ(alphabet_size_v<TypeParam>,
        static_cast<size_t>(std::numeric_limits<TypeParam>::max()) + 1 - std::numeric_limits<TypeParam>::lowest());
}
