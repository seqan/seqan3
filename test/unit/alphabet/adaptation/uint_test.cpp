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

// uint32_t, too slow
using fast_uint_types = ::testing::Types<uint8_t, uint16_t/*, uint32_t*/>;

INSTANTIATE_TYPED_TEST_CASE_P(uint_adaptation, alphabet, fast_uint_types);
INSTANTIATE_TYPED_TEST_CASE_P(uint_adaptation, alphabet_constexpr, fast_uint_types);

template <typename T>
class uint_adaptation : public ::testing::Test
{};

using uint_types = ::testing::Types<uint8_t, uint16_t, uint32_t>;

TYPED_TEST_CASE(uint_adaptation, uint_types);

TYPED_TEST(uint_adaptation, concept_check)
{
    EXPECT_TRUE(uint_adaptation_concept<TypeParam   >);
    EXPECT_TRUE(uint_adaptation_concept<TypeParam & >);
    EXPECT_TRUE(uint_adaptation_concept<TypeParam &&>);
}

TYPED_TEST(uint_adaptation, underlying_rank_t)
{
    EXPECT_TRUE((std::is_same_v<underlying_rank_t<TypeParam   >, TypeParam>));
    EXPECT_TRUE((std::is_same_v<underlying_rank_t<TypeParam & >, TypeParam>));
    EXPECT_TRUE((std::is_same_v<underlying_rank_t<TypeParam &&>, TypeParam>));
}

TYPED_TEST(uint_adaptation, to_rank)
{
    TypeParam l{65};
    EXPECT_TRUE((std::is_same_v<decltype(to_rank(l)),             underlying_rank_t<TypeParam>>));
    EXPECT_TRUE((std::is_same_v<decltype(to_rank(TypeParam{65})), underlying_rank_t<TypeParam>>));
    EXPECT_EQ(to_rank(TypeParam{65}), l);
}

TYPED_TEST(uint_adaptation, assign_rank)
{
    TypeParam l{65};
    EXPECT_TRUE((std::is_same_v<decltype(assign_rank(l,             65)), underlying_rank_t<TypeParam> &>));
    EXPECT_TRUE((std::is_same_v<decltype(assign_rank(TypeParam{65}, 65)), underlying_rank_t<TypeParam>  >));
    EXPECT_EQ((assign_rank(TypeParam{65}, 65)), l);
    EXPECT_EQ((assign_rank(l,             67)), TypeParam{67});
}

TYPED_TEST(uint_adaptation, underlying_char_t)
{
    EXPECT_TRUE((std::is_integral_v<underlying_rank_t<TypeParam>>));
    EXPECT_GE(sizeof(underlying_rank_t<TypeParam>), sizeof(TypeParam));
}

TYPED_TEST(uint_adaptation, to_char)
{
    TypeParam l{65};
    EXPECT_TRUE((std::is_same_v<decltype(to_char(l)),             underlying_char_t<TypeParam>>));
    EXPECT_TRUE((std::is_same_v<decltype(to_char(TypeParam{65})), underlying_char_t<TypeParam>>));
    if constexpr (std::is_unsigned_v<TypeParam>)
    {
        unsigned char cmp{'A'};
        EXPECT_EQ(to_rank(TypeParam{65}), cmp);
    }
    else
    {
        EXPECT_EQ(to_rank(TypeParam{65}), 'A');
    }
}

TYPED_TEST(uint_adaptation, assign_char)
{
    TypeParam l{65};
    EXPECT_TRUE((std::is_same_v<decltype(assign_char(l,              'A')), underlying_rank_t<TypeParam> &>));
    EXPECT_TRUE((std::is_same_v<decltype(assign_char(TypeParam{'A'}, 'A')), underlying_rank_t<TypeParam>  >));
    EXPECT_EQ((assign_char(TypeParam{67}, 'A')), l);
    EXPECT_EQ((assign_char(l,             'C')), TypeParam{67});
}

TYPED_TEST(uint_adaptation, assign_char_strict)
{
    TypeParam l{65};
    EXPECT_TRUE((std::is_same_v<decltype(assign_char_strict(l,              'A')), underlying_rank_t<TypeParam> &>));
    EXPECT_TRUE((std::is_same_v<decltype(assign_char_strict(TypeParam{'A'}, 'A')), underlying_rank_t<TypeParam>  >));
    EXPECT_EQ((assign_char_strict(TypeParam{67}, 'A')), l);
    EXPECT_EQ((assign_char_strict(l,             'C')), TypeParam{67});
}

TYPED_TEST(uint_adaptation, alphabet_size_v)
{
    EXPECT_EQ(alphabet_size_v<TypeParam>,
        static_cast<uint64_t>(std::numeric_limits<TypeParam>::max()) + 1 - std::numeric_limits<TypeParam>::lowest());
}
