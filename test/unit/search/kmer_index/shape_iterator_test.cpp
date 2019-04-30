// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/kmer_index/shape_iterator.hpp>

using namespace seqan3;

// ============================================================================
//  shape
// ============================================================================

TEST(shape, ctr)
{
    EXPECT_TRUE((std::is_default_constructible_v<shape>));
    EXPECT_TRUE((std::is_nothrow_default_constructible_v<shape>));
    EXPECT_TRUE((std::is_copy_constructible_v<shape>));
    EXPECT_TRUE((std::is_trivially_copy_constructible_v<shape>));
    EXPECT_TRUE((std::is_nothrow_copy_constructible_v<shape>));
    EXPECT_TRUE((std::is_move_constructible_v<shape>));
    EXPECT_TRUE((std::is_trivially_move_constructible_v<shape>));
    EXPECT_TRUE((std::is_nothrow_move_constructible_v<shape>));
    EXPECT_TRUE((std::is_copy_assignable_v<shape>));
    EXPECT_TRUE((std::is_trivially_copy_assignable_v<shape>));
    EXPECT_TRUE((std::is_nothrow_copy_assignable_v<shape>));
    EXPECT_TRUE((std::is_move_assignable_v<shape>));
    EXPECT_TRUE((std::is_trivially_move_assignable_v<shape>));
    EXPECT_TRUE((std::is_nothrow_move_assignable_v<shape>));

    shape s0(0, 0, 0);
    std::array<bool, 32> a0{false, false, false};
    EXPECT_EQ(s0, a0);

    shape s1(1, 0, 1);
    std::array<bool, 32> a1{true, false, true};
    EXPECT_EQ(s1, a1);

    shape s2(3);
    std::array<bool, 32> a2{true, true, true};
    EXPECT_EQ(s2, a2);
}

TEST(shape, size)
{
    EXPECT_EQ(std::ranges::size(shape(1)), 1u);
    EXPECT_EQ(std::ranges::size(shape(32)), 32u);
    EXPECT_EQ(std::ranges::size(shape(0, 0)), 2u);
    EXPECT_EQ(std::ranges::size(shape(1, 0, 1, 0, 1, 1)), 6u);
}

// ============================================================================
//  iterator
// ============================================================================


TEST(shape_iterator, ctr)
{
    auto text1{"C"_dna4};
    shape_iterator it(std::ranges::begin(text1), shape(1));
    EXPECT_EQ(*it, 1u);

    auto text2{"ACGT"_dna4};
    shape_iterator it2(std::ranges::begin(text2), shape(4));
    EXPECT_EQ(*it2, 27u);

    shape_iterator it3(std::ranges::begin(text2), shape(0, 0, 0, 1));
    EXPECT_EQ(*it3, 3u);
}

TEST(shape_iterator, increment)
{
    std::array<size_t, 3> expected1{1, 6, 11};
    auto text1{"ACGT"_dna4};
    shape_iterator it1(std::ranges::begin(text1), shape(2));

    EXPECT_EQ(*it1, expected1[0]);
    EXPECT_EQ(*(++it1), expected1[1]);
    EXPECT_EQ(*(++it1), expected1[2]);

    std::array<size_t, 3> expected2{2, 19, 32};
    auto text2{"ACGTA"_dna4};
    shape_iterator it2(std::ranges::begin(text2), shape(1, 0, 1));

    EXPECT_EQ(*it2, expected2[0]);
    EXPECT_EQ(*(++it2), expected2[1]);
    EXPECT_EQ(*(++it2), expected2[2]);
}

TEST(shape_iterator, random_access)
{
    std::array<size_t, 3> expected{2, 19, 32};
    auto text{"ACGTA"_dna4};
    shape_iterator it(std::ranges::begin(text), shape(1, 0, 1));

    EXPECT_EQ(*it[2], expected[2]);
    EXPECT_EQ(*it[0], expected[0]);
    EXPECT_EQ(*it[1], expected[1]);
}

TEST(shape_iterator, comparison)
{
    auto text{"ACGT"_dna4};
    shape_iterator it1(std::ranges::begin(text), shape(2));
    shape_iterator it2(std::next(std::ranges::begin(text), 2), shape(2));

    EXPECT_NE(it1, it2);
    EXPECT_NE(it1, std::ranges::begin(text));
    EXPECT_EQ(++it1, it2);
    EXPECT_EQ(it1, std::next(std::ranges::begin(text), 2));
}

TEST(shape_iterator, sentinel)
{
    auto text{"ACGT"_dna4};
    size_t i{0};

    for (shape_iterator it(std::ranges::begin(text), shape(2)); it != std::ranges::end(text); ++it)
        ++i;

    EXPECT_EQ(i, 3u);
}
