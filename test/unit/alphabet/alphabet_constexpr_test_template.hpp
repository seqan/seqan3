// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/concept.hpp>

template <typename T>
using alphabet_constexpr = ::testing::Test;

TYPED_TEST_SUITE_P(alphabet_constexpr);

TYPED_TEST_P(alphabet_constexpr, concept_check)
{
    EXPECT_TRUE(seqan3::detail::constexpr_alphabet<TypeParam   >);
    EXPECT_TRUE(seqan3::detail::constexpr_alphabet<TypeParam & >);

    EXPECT_TRUE(seqan3::detail::constexpr_alphabet<TypeParam const   >);
    EXPECT_TRUE(seqan3::detail::constexpr_alphabet<TypeParam const & >);

    EXPECT_TRUE(seqan3::detail::writable_constexpr_alphabet<TypeParam   >);
    EXPECT_TRUE(seqan3::detail::writable_constexpr_alphabet<TypeParam & >);

    EXPECT_FALSE(seqan3::detail::writable_constexpr_alphabet<TypeParam const   >);
    EXPECT_FALSE(seqan3::detail::writable_constexpr_alphabet<TypeParam const & >);
}

TYPED_TEST_P(alphabet_constexpr, global_assign_char)
{
    [[maybe_unused]] constexpr TypeParam t0{seqan3::assign_char_to('A', TypeParam{})};
}

TYPED_TEST_P(alphabet_constexpr, global_to_char)
{
    constexpr TypeParam t0{TypeParam{}};
    [[maybe_unused]] constexpr seqan3::alphabet_char_t<TypeParam> c = seqan3::to_char(t0);
}

REGISTER_TYPED_TEST_SUITE_P(alphabet_constexpr,
                            concept_check,
                            global_assign_char,
                            global_to_char);
