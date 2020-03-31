// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/adaptation/all.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

// uint32_t, too slow
using fast_uint_types = ::testing::Types<uint8_t, uint16_t/*, uint32_t*/>;

INSTANTIATE_TYPED_TEST_SUITE_P(uint_adaptation, alphabet_, fast_uint_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(uint_adaptation, semi_alphabet_test, fast_uint_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(uint_adaptation, alphabet_constexpr, fast_uint_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(uint_adaptation, semi_alphabet_constexpr, fast_uint_types, );

template <typename T>
using uint_adaptation = ::testing::Test;

using uint_types = ::testing::Types<uint8_t, uint16_t, uint32_t>;

TYPED_TEST_SUITE(uint_adaptation, uint_types, );

TYPED_TEST(uint_adaptation, type_properties)
{
    EXPECT_TRUE((std::is_trivially_copyable_v<TypeParam>));
    EXPECT_TRUE((std::is_trivially_default_constructible_v<TypeParam>));
    EXPECT_TRUE((std::is_trivial_v<TypeParam>));
}

TYPED_TEST(uint_adaptation, alphabet_rank_t)
{
    EXPECT_TRUE((std::is_same_v<seqan3::alphabet_rank_t<TypeParam>, TypeParam>));
    EXPECT_TRUE((std::is_same_v<seqan3::alphabet_rank_t<TypeParam &>, TypeParam>));
    EXPECT_TRUE((std::is_same_v<seqan3::alphabet_rank_t<TypeParam &&>, TypeParam>));
}

TYPED_TEST(uint_adaptation, to_rank)
{
    TypeParam l{65};
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::to_rank(l)),
                                seqan3::alphabet_rank_t<TypeParam>>));
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::to_rank(TypeParam{65})),
                                seqan3::alphabet_rank_t<TypeParam>>));
    EXPECT_EQ(seqan3::to_rank(TypeParam{65}), l);
}

TYPED_TEST(uint_adaptation, assign_rank)
{
    TypeParam l{65};
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::assign_rank_to(65, l)),
                                seqan3::alphabet_rank_t<TypeParam> &>));
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::assign_rank_to(65, TypeParam{65})),
                                seqan3::alphabet_rank_t<TypeParam>>));
    EXPECT_EQ((seqan3::assign_rank_to(65, TypeParam{65})), l);
    EXPECT_EQ((seqan3::assign_rank_to(67, l)), TypeParam{67});
}

TYPED_TEST(uint_adaptation, alphabet_char_t)
{
    EXPECT_TRUE((std::is_integral_v<seqan3::alphabet_rank_t<TypeParam>>));
    EXPECT_GE(sizeof(seqan3::alphabet_rank_t<TypeParam>), sizeof(TypeParam));
}

TYPED_TEST(uint_adaptation, to_char)
{
    TypeParam l{65};
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::to_char(l)),
                                seqan3::alphabet_char_t<TypeParam>>));
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::to_char(TypeParam{65})),
                                seqan3::alphabet_char_t<TypeParam>>));
    if constexpr (std::is_unsigned_v<TypeParam>)
    {
        unsigned char cmp{'A'};
        EXPECT_EQ(seqan3::to_rank(TypeParam{65}), cmp);
    }
    else
    {
        EXPECT_EQ(seqan3::to_rank(TypeParam{65}), 'A');
    }
}

TYPED_TEST(uint_adaptation, assign_char)
{
    TypeParam l{65};
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::assign_char_to('A', l)),
                                seqan3::alphabet_rank_t<TypeParam> &>));
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::assign_char_to('A', TypeParam{'A'})),
                                seqan3::alphabet_rank_t<TypeParam>>));
    EXPECT_EQ((seqan3::assign_char_to('A', TypeParam{67})), l);
    EXPECT_EQ((seqan3::assign_char_to('C', l)), TypeParam{67});
}

TYPED_TEST(uint_adaptation, assign_char_strictly_to)
{
    TypeParam l{65};
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::assign_char_strictly_to('A', l)),
                                seqan3::alphabet_rank_t<TypeParam> &>));
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::assign_char_strictly_to('A', TypeParam{'A'})),
                                seqan3::alphabet_rank_t<TypeParam>>));
    EXPECT_EQ((seqan3::assign_char_strictly_to('A', TypeParam{67})), l);
    EXPECT_EQ((seqan3::assign_char_strictly_to('C', l)), TypeParam{67});
}

TYPED_TEST(uint_adaptation, alphabet_size)
{
    EXPECT_EQ(seqan3::alphabet_size<TypeParam>,
        static_cast<uint64_t>(std::numeric_limits<TypeParam>::max()) + 1 - std::numeric_limits<TypeParam>::lowest());
}
