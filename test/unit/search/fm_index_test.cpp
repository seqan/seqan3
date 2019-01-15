// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/test/tmp_filename.hpp>

#include <gtest/gtest.h>

using namespace seqan3;

// TODO: EXPECT_EQ is not supported by sdsl

template <typename T>
class fm_index_test : public ::testing::Test
{};

using fm_index_types = ::testing::Types<fm_index<std::vector<dna4>>,
                                        bi_fm_index<std::vector<dna4>>,
                                        bi_fm_index<std::vector<aa27>>,
                                        bi_fm_index<std::vector<char>>>;

TYPED_TEST_CASE(fm_index_test, fm_index_types);

TYPED_TEST(fm_index_test, ctr)
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

TYPED_TEST(fm_index_test, swap)
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

TYPED_TEST(fm_index_test, size)
{
    TypeParam fm;
    EXPECT_TRUE(fm.empty());

    typename TypeParam::text_type test(8);
    fm.construct(test);
    EXPECT_EQ(fm.size(), 9u); // including a sentinel character
}

TYPED_TEST(fm_index_test, serialization)
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

TEST(fm_index_test, concepts)
{
    EXPECT_TRUE(FmIndex<fm_index<std::vector<dna4>>>);
    EXPECT_TRUE(FmIndex<fm_index<std::vector<dna5>>>);
    EXPECT_TRUE(FmIndexTraits<fm_index_default_traits>);

    EXPECT_TRUE(BiFmIndex<bi_fm_index<std::vector<dna4>>>);
    EXPECT_TRUE(BiFmIndex<bi_fm_index<std::vector<dna5>>>);
    EXPECT_TRUE(BiFmIndexTraits<bi_fm_index_default_traits>);
}
