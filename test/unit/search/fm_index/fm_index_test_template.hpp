// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/fm_index/concept.hpp>
#include <seqan3/test/cereal.hpp>

template <typename T>
class fm_index_test : public ::testing::Test
{};

TYPED_TEST_SUITE_P(fm_index_test);

TYPED_TEST_P(fm_index_test, ctr)
{
    using index_t = typename TypeParam::first_type;
    using text_t = typename TypeParam::second_type;

    text_t text(10); // initialized with smallest char

    // construction from text
    index_t fm0{text};

    // copy construction
    index_t fm1{fm0};
    EXPECT_EQ(fm0, fm1);

    // copy assignment
    index_t fm2 = fm0;
    EXPECT_EQ(fm0, fm2);

    // move construction
    index_t fm3{std::move(fm1)};
    EXPECT_EQ(fm0, fm3);

    // move assigment
    index_t fm4 = std::move(fm2);
    EXPECT_EQ(fm0, fm4);

    // container contructor
    index_t fm5{text};
    EXPECT_EQ(fm0, fm5);
}

TYPED_TEST_P(fm_index_test, swap)
{
    using index_t = typename TypeParam::first_type;
    using text_t = typename TypeParam::second_type;

    text_t textA(10);
    text_t textB(20);

    index_t fm0{textA};
    index_t fm1{textB};
    index_t fm2{fm0};
    index_t fm3{fm1};

    EXPECT_EQ(fm0, fm2);
    EXPECT_EQ(fm1, fm3);
    EXPECT_NE(fm0, fm1);

    std::swap(fm1, fm2);

    EXPECT_EQ(fm0, fm1);
    EXPECT_EQ(fm2, fm3);
    EXPECT_NE(fm0, fm2);
}

TYPED_TEST_P(fm_index_test, size)
{
    using index_t = typename TypeParam::first_type;
    using text_t = typename TypeParam::second_type;

    index_t fm;
    EXPECT_TRUE(fm.empty());

    text_t text(8);
    fm = index_t{text};
    EXPECT_EQ(fm.size(), 9u); // including a sentinel character
}

TYPED_TEST_P(fm_index_test, empty_text)
{
    using index_t = typename TypeParam::first_type;
    using text_t = typename TypeParam::second_type;

    text_t text{};
    EXPECT_THROW(index_t index{text}, std::invalid_argument);
}

TYPED_TEST_P(fm_index_test, serialisation)
{
    using index_t = typename TypeParam::first_type;
    using text_t = typename TypeParam::second_type;

    text_t text(10);

    index_t fm{text};
    seqan3::test::do_serialisation(fm);
}

REGISTER_TYPED_TEST_SUITE_P(fm_index_test, ctr, swap, size, empty_text, serialisation);
