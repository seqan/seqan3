/// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/test/cereal.hpp>

template <typename T>
class fm_index_collection_test : public ::testing::Test
{};

TYPED_TEST_SUITE_P(fm_index_collection_test);

TYPED_TEST_P(fm_index_collection_test, ctr)
{
    using index_t = typename TypeParam::first_type;
    using text_t = typename TypeParam::second_type;
    using inner_text_type = std::ranges::range_value_t<text_t>;

    text_t text{inner_text_type(10), inner_text_type(10)}; // initialized with smallest char

    // default/zero construction and construct from text
    index_t fm0{text}, fm2, fm4;

    // copy construction
    index_t fm1{fm0};
    EXPECT_EQ(fm0, fm1);
    // Make sure rank and select support pointer are correct by using them.
    auto it0 = fm0.cursor();
    it0.extend_right(inner_text_type(5));
    auto it1 = fm1.cursor();
    it1.extend_right(inner_text_type(5));
    EXPECT_EQ(it0.locate(), it1.locate());

    // copy assignment
    fm2 = fm0;
    EXPECT_EQ(fm0, fm2);
    // Make sure rank and select support pointer are correct by using them.
    auto it2 = fm2.cursor();
    it2.extend_right(inner_text_type(5));
    EXPECT_EQ(it0.locate(), it2.locate());

    // move construction
    index_t fm3{std::move(fm1)};
    EXPECT_EQ(fm0, fm3);
    // Make sure rank and select support pointer are correct by using them.
    auto it3 = fm3.cursor();
    it3.extend_right(inner_text_type(5));
    EXPECT_EQ(it0.locate(), it3.locate());

    // move assigment
    fm4 = std::move(fm2);
    EXPECT_EQ(fm0, fm4);
    // Make sure rank and select support pointer are correct by using them.
    auto it4 = fm4.cursor();
    it4.extend_right(inner_text_type(5));
    EXPECT_EQ(it0.locate(), it4.locate());

    // container contructor
    index_t fm5{text};
    EXPECT_EQ(fm0, fm5);
}

TYPED_TEST_P(fm_index_collection_test, swap)
{
    using index_t = typename TypeParam::first_type;
    using text_t = typename TypeParam::second_type;
    using inner_text_type = std::ranges::range_value_t<text_t>;

    text_t textA{inner_text_type(10), inner_text_type(10)};
    text_t textB{inner_text_type(20), inner_text_type(20)};

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

    std::swap(fm0, fm1);
    // Make sure rank and select support pointer are correct by using them.
    auto it0 = fm0.cursor();
    it0.extend_right(inner_text_type(5));
    auto it1 = fm1.cursor();
    it1.extend_right(inner_text_type(5));
    EXPECT_EQ(it0.locate(), it1.locate());
}

TYPED_TEST_P(fm_index_collection_test, size)
{
    using index_t = typename TypeParam::first_type;
    using text_t = typename TypeParam::second_type;
    using inner_text_type = std::ranges::range_value_t<text_t>;

    index_t fm;
    EXPECT_TRUE(fm.empty());

    text_t text{inner_text_type(4), inner_text_type(4)};
    fm = index_t{text};
    EXPECT_EQ(fm.size(), 10u); // including a sentinel character and a delimiter
}

TYPED_TEST_P(fm_index_collection_test, concept_check)
{
    using index_t = typename TypeParam::first_type;

    EXPECT_TRUE((seqan3::fm_index_specialisation<index_t>));
    if constexpr (std::same_as<index_t, seqan3::bi_fm_index<typename index_t::alphabet_type,
                                                            seqan3::text_layout::collection>>)
    {
        EXPECT_TRUE(seqan3::bi_fm_index_specialisation<index_t>);
    }
}

TYPED_TEST_P(fm_index_collection_test, empty_text)
{
    using index_t = typename TypeParam::first_type;
    using text_t = typename TypeParam::second_type;

    {
        text_t text{};
        EXPECT_THROW(index_t index{text}, std::invalid_argument);
    }
    {
        text_t text(2);
        EXPECT_THROW(index_t index{text}, std::invalid_argument);
    }
}

TYPED_TEST_P(fm_index_collection_test, serialisation)
{
    using index_t = typename TypeParam::first_type;
    using text_t = typename TypeParam::second_type;
    using inner_text_type = std::ranges::range_value_t<text_t>;

    text_t text{inner_text_type(4), inner_text_type(12)};

    index_t fm{text};
    seqan3::test::do_serialisation(fm);
}

REGISTER_TYPED_TEST_SUITE_P(fm_index_collection_test, ctr, swap, size, serialisation, concept_check, empty_text);
