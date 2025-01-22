// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/concept.hpp>

template <typename t>
using alphabet_constexpr = ::testing::Test;

TYPED_TEST_SUITE_P(alphabet_constexpr);

TYPED_TEST_P(alphabet_constexpr, concept_check)
{
    EXPECT_TRUE(seqan3::detail::constexpr_alphabet<TypeParam>);
    EXPECT_TRUE(seqan3::detail::constexpr_alphabet<TypeParam &>);

    EXPECT_TRUE(seqan3::detail::constexpr_alphabet<TypeParam const>);
    EXPECT_TRUE(seqan3::detail::constexpr_alphabet<TypeParam const &>);

    EXPECT_TRUE(seqan3::detail::writable_constexpr_alphabet<TypeParam>);
    EXPECT_TRUE(seqan3::detail::writable_constexpr_alphabet<TypeParam &>);

    EXPECT_FALSE(seqan3::detail::writable_constexpr_alphabet<TypeParam const>);
    EXPECT_FALSE(seqan3::detail::writable_constexpr_alphabet<TypeParam const &>);
}

TYPED_TEST_P(alphabet_constexpr, assign_char)
{
    [[maybe_unused]] constexpr TypeParam t0{seqan3::assign_char_to('A', TypeParam{})};
}

TYPED_TEST_P(alphabet_constexpr, to_char)
{
    constexpr TypeParam t0{TypeParam{}};
    [[maybe_unused]] constexpr seqan3::alphabet_char_t<TypeParam> c = seqan3::to_char(t0);
}

REGISTER_TYPED_TEST_SUITE_P(alphabet_constexpr, concept_check, assign_char, to_char);
