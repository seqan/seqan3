/// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/test/cereal.hpp>

using namespace seqan3;

template <typename T>
class fm_index_collection_test : public ::testing::Test
{};

TYPED_TEST_CASE_P(fm_index_collection_test);

TYPED_TEST_P(fm_index_collection_test, ctr)
{
    using inner_text_type = value_type_t<typename TypeParam::text_type>;
    typename TypeParam::text_type text{inner_text_type(10), inner_text_type(10)}; // initialized with smallest char

    // default/zero construction
    TypeParam fm0;
    fm0.construct(text);

    // copy construction
    TypeParam fm1{fm0};
    EXPECT_EQ(fm0, fm1);

    // copy assignment
    TypeParam fm2 = fm0;
    EXPECT_EQ(fm0, fm2);

    // move construction
    TypeParam fm3{std::move(fm1)};
    EXPECT_EQ(fm0, fm3);

    // move assigment
    TypeParam fm4 = std::move(fm2);
    EXPECT_EQ(fm0, fm4);

    // container contructor
    TypeParam fm5{text};
    EXPECT_EQ(fm0, fm5);
}

TYPED_TEST_P(fm_index_collection_test, swap)
{
    using inner_text_type = value_type_t<typename TypeParam::text_type>;
    typename TypeParam::text_type textA{inner_text_type(10), inner_text_type(10)};
    typename TypeParam::text_type textB{inner_text_type(20), inner_text_type(20)};

    TypeParam fm0{textA};
    TypeParam fm1{textB};
    TypeParam fm2{fm0};
    TypeParam fm3{fm1};

    EXPECT_EQ(fm0, fm2);
    EXPECT_EQ(fm1, fm3);
    EXPECT_NE(fm0, fm1);

    std::swap(fm1, fm2);

    EXPECT_EQ(fm0, fm1);
    EXPECT_EQ(fm2, fm3);
    EXPECT_NE(fm0, fm2);
}

TYPED_TEST_P(fm_index_collection_test, size)
{
    TypeParam fm;
    EXPECT_TRUE(fm.empty());

    using inner_text_type = value_type_t<typename TypeParam::text_type>;
    typename TypeParam::text_type text{inner_text_type(4), inner_text_type(4)};
    fm.construct(text);
    EXPECT_EQ(fm.size(), 10u); // including a sentinel character and a delimiter
}

TYPED_TEST_P(fm_index_collection_test, concept_check)
{
    EXPECT_TRUE(FmIndex<TypeParam>);
    if constexpr (detail::is_type_specialisation_of_v<TypeParam, bi_fm_index>)
    {
        EXPECT_TRUE(BiFmIndex<TypeParam>);
    }
}

TYPED_TEST_P(fm_index_collection_test, empty_text)
{
    {
        typename TypeParam::text_type text{};
        EXPECT_THROW(TypeParam index{text}, std::invalid_argument);
    }
    {
        typename TypeParam::text_type text(2);
        EXPECT_THROW(TypeParam index{text}, std::invalid_argument);
    }
}

TYPED_TEST_P(fm_index_collection_test, serialisation)
{
    using inner_text_type = value_type_t<typename TypeParam::text_type>;
    typename TypeParam::text_type text{inner_text_type(4), inner_text_type(12)};

    TypeParam fm{text};
    test::do_serialisation(fm);
}

REGISTER_TYPED_TEST_CASE_P(fm_index_collection_test, ctr, swap, size, serialisation, concept_check, empty_text);
