/// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <gtest/gtest.h>

using namespace seqan3;

// TODO: EXPECT_EQ is not supported by sdsl

template <typename T>
class fm_index_test : public ::testing::Test
{};

TYPED_TEST_CASE_P(fm_index_test);

TYPED_TEST_P(fm_index_test, ctr)
{
    typename TypeParam::text_type text(10); // initialized with smallest char

    // default/zero construction
    TypeParam fm0;
    fm0.construct(text);

    // copy construction
    TypeParam fm1{fm0};
    EXPECT_EQ(fm0.size(), fm1.size());

    // copy assignment
    TypeParam fm2 = fm0;
    EXPECT_EQ(fm0.size(), fm2.size());

    // move construction
    TypeParam fm3{std::move(fm0)};
    EXPECT_EQ(fm0.size(), fm3.size());

    // move assigment
    TypeParam fm4 = std::move(fm0);
    EXPECT_EQ(fm0.size(), fm4.size());

    // container contructor
    TypeParam fm5{text};
    EXPECT_EQ(fm0.size(), fm5.size());
}

TYPED_TEST_P(fm_index_test, swap)
{
    typename TypeParam::text_type textA(10);
    typename TypeParam::text_type textB(20);

    TypeParam fm0{textA};
    TypeParam fm1{textB};
    TypeParam fm2{fm0};
    TypeParam fm3{fm1};

    EXPECT_EQ(fm0.size(), fm2.size());
    EXPECT_EQ(fm1.size(), fm3.size());
    EXPECT_NE(fm0.size(), fm1.size());

    std::swap(fm1, fm2);

    EXPECT_EQ(fm0.size(), fm1.size());
    EXPECT_EQ(fm2.size(), fm3.size());
    EXPECT_NE(fm0.size(), fm2.size());
}

TYPED_TEST_P(fm_index_test, size)
{
    TypeParam fm;
    EXPECT_TRUE(fm.empty());

    typename TypeParam::text_type test(8);
    fm.construct(test);
    EXPECT_EQ(fm.size(), 9u); // including a sentinel character
}

TYPED_TEST_P(fm_index_test, serialization)
{
    typename TypeParam::text_type text(8);
    TypeParam fm0{text};

    test::tmp_filename filename{"fm_index"};
    auto const & path = filename.get_path();

    EXPECT_TRUE(fm0.store(path));

    TypeParam fm1{};
    EXPECT_TRUE(fm1.load(path));

    EXPECT_EQ(fm1.size(), 9u);
}

TYPED_TEST_P(fm_index_test, concept_check)
{
    EXPECT_TRUE(FmIndex<TypeParam>);
    if constexpr (detail::is_type_specialisation_of_v<TypeParam, bi_fm_index>)
    {
        EXPECT_TRUE(BiFmIndex<TypeParam>);
    }
}

TYPED_TEST_P(fm_index_test, empty_text)
{
    typename TypeParam::text_type text{};
    EXPECT_THROW(TypeParam index{text}, std::invalid_argument);
}

REGISTER_TYPED_TEST_CASE_P(fm_index_test, ctr, swap, size, serialization, concept_check, empty_text);
