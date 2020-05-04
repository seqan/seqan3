/// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>
#include <type_traits>

#include <seqan3/range/views/slice.hpp>
#include <seqan3/search/fm_index/bi_fm_index_cursor.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/test/expect_range_eq.hpp>

#include "../helper.hpp"

template <typename T>
struct bi_fm_index_cursor_test;

TYPED_TEST_SUITE_P(bi_fm_index_cursor_test);

TYPED_TEST_P(bi_fm_index_cursor_test, cursor)
{
    typename TypeParam::index_type bi_fm{this->text1};  // "AACGATCGGA"

    seqan3::fm_index fm_fwd{this->text1};
    seqan3::fm_index fm_rev{this->rev_text2};

    TypeParam bi_it = bi_fm.cursor();
    EXPECT_EQ(seqan3::uniquify(bi_it.locate()), seqan3::uniquify(bi_fm.fwd_begin().locate()));
    EXPECT_EQ(seqan3::uniquify(bi_it.locate()), seqan3::uniquify(bi_fm.rev_begin().locate()));
}

TYPED_TEST_P(bi_fm_index_cursor_test, extend)
{
    typename TypeParam::index_type bi_fm{seqan3::views::slice(this->text, 0, 11)};  // "ACGGTAGGACG"

    auto it = bi_fm.cursor();
    EXPECT_TRUE(it.extend_right()); // "A"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{0, 5, 8}));
    EXPECT_TRUE(it.extend_left()); // "GA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{7}));
    EXPECT_TRUE(it.extend_right()); // "GAC"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{7}));
    EXPECT_TRUE(it.extend_right()); // "GACG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{7}));
    EXPECT_FALSE(it.extend_right()); // "GACG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{7}));
    EXPECT_TRUE(it.extend_left()); // "GGACG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{6}));
}

TYPED_TEST_P(bi_fm_index_cursor_test, extend_char)
{
    typename TypeParam::index_type bi_fm{seqan3::views::slice(this->text, 0, 11)};  // "ACGGTAGGACG"

    auto it = bi_fm.cursor();
    EXPECT_TRUE(it.extend_left(this->text[2])); // "G"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{2, 3, 6, 7, 10}));
    EXPECT_TRUE(it.extend_left(this->text[1])); // "CG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{1, 9}));
    EXPECT_FALSE(it.extend_left(this->text[1])); // "CG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{1, 9}));
    EXPECT_FALSE(it.extend_left(this->text[2])); // "CG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{1, 9}));
    EXPECT_FALSE(it.extend_right(this->text[4])); // "CG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{1, 9}));
    EXPECT_TRUE(it.extend_right(this->text[2])); // "CGG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{1}));
    EXPECT_TRUE(it.extend_right(this->text[4])); // "CGGT"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{1}));
    EXPECT_TRUE(it.extend_right(this->text[0])); // "CGGTA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{1}));
    EXPECT_TRUE(it.extend_left(this->text[0])); // "ACGGTA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{0}));
    EXPECT_FALSE(it.extend_left(this->text[0])); // "ACGGTA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{0}));
}

TYPED_TEST_P(bi_fm_index_cursor_test, extend_range)
{
    typename TypeParam::index_type bi_fm{seqan3::views::slice(this->text, 0, 11)};  // "ACGGTAGGACG"

    auto it = bi_fm.cursor();
    EXPECT_FALSE(it.extend_left(this->pattern1)); // ""
    // sentinel position included
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
    EXPECT_TRUE(it.extend_left(seqan3::views::slice(this->text, 1, 3))); // "CG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{1, 9}));
    EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text, 3, 6))); // "CGGTA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{1}));
    EXPECT_FALSE(it.extend_left(seqan3::views::slice(this->text, 4, 6))); // "CGGTA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{1}));
    EXPECT_TRUE(it.extend_left(seqan3::views::slice(this->text, 0, 1))); // "ACGGTA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{0}));
}

TYPED_TEST_P(bi_fm_index_cursor_test, extend_and_cycle)
{
    typename TypeParam::index_type bi_fm{seqan3::views::slice(this->text, 0, 11)};  // "ACGGTAGGACG"

    auto it = bi_fm.cursor();
    EXPECT_TRUE(it.extend_right()); // "A"
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_front(), "");
#endif
    EXPECT_TRUE(it.extend_left()); // "GA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{7}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_back(), "");
#endif
    EXPECT_TRUE(it.cycle_front()); // "TA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{4}));
    EXPECT_FALSE(it.cycle_front()); // "TA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{4}));
}

TYPED_TEST_P(bi_fm_index_cursor_test, extend_range_and_cycle)
{
    typename TypeParam::index_type bi_fm{seqan3::views::slice(this->text, 0, 14)};  // "ACGGTAGGACGTAG"

    auto it = bi_fm.cursor();
    EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text, 0, 2))); // "AC"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{0, 8}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_front(), "");
#endif
    EXPECT_TRUE(it.cycle_back()); // "AG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{5, 12}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_front(), "");
#endif
    EXPECT_FALSE(it.extend_left(this->pattern2)); // "AG" ("TT")
    EXPECT_TRUE(it.extend_left(seqan3::views::slice(this->text, 9, 12))); // "CGTAG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{9}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_back(), "");
#endif
    EXPECT_TRUE(it.cycle_front()); // "GGTAG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{2}));
}

TYPED_TEST_P(bi_fm_index_cursor_test, to_fwd_cursor)
{
    typename TypeParam::index_type bi_fm{this->text};   // "ACGGTAGGACGTAGC"

    {
        auto it = bi_fm.cursor();
        EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text, 10, 15))); // "GTAGC"
        EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{10}));

        auto fwd_it = it.to_fwd_cursor();
        EXPECT_TRUE(fwd_it.cycle_back()); // "GTAGG"
        EXPECT_EQ(seqan3::uniquify(fwd_it.locate()), (std::vector<uint64_t>{3}));
        EXPECT_RANGE_EQ(fwd_it.path_label(this->text), seqan3::views::slice(this->text, 3, 8)); // "GTAGG"
        EXPECT_FALSE(fwd_it.cycle_back());
    }

    {
        auto it = bi_fm.cursor();
        EXPECT_TRUE(it.extend_left(seqan3::views::slice(this->text, 3, 7))); // "GTAG"
        EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{3, 10}));

        auto fwd_it = it.to_fwd_cursor();
    #ifndef NDEBUG
        EXPECT_DEATH(fwd_it.cycle_back(), "");
    #endif
        EXPECT_TRUE(fwd_it.extend_right());
        EXPECT_EQ(seqan3::uniquify(fwd_it.locate()), (std::vector<uint64_t>{10}));
        EXPECT_RANGE_EQ(fwd_it.path_label(this->text), seqan3::views::slice(this->text, 10, 15)); // "GTAGC"
        EXPECT_TRUE(fwd_it.cycle_back());
        EXPECT_EQ(seqan3::uniquify(fwd_it.locate()), (std::vector<uint64_t>{3}));
        EXPECT_RANGE_EQ(fwd_it.path_label(this->text), seqan3::views::slice(this->text, 3, 8)); // "GTAGG"
    }
}

TYPED_TEST_P(bi_fm_index_cursor_test, to_rev_cursor)
{

    typename TypeParam::index_type bi_fm{this->text};   // "ACGGTAGGACGTAGC"

    {
        auto it = bi_fm.cursor();
        EXPECT_TRUE(it.extend_left(seqan3::views::slice(this->text, 9, 14))); // "CGTAG"
        EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{9}));

        auto rev_it = it.to_rev_cursor(); // text "CGATGCAGGATGGCA"
        EXPECT_EQ(seqan3::uniquify(rev_it.locate()), (std::vector<uint64_t>{1}));
        EXPECT_RANGE_EQ(rev_it.path_label(this->rev_text1), this->pattern3);    // "GATGC"
        EXPECT_TRUE(rev_it.cycle_back()); // "GATGG"
        EXPECT_EQ(seqan3::uniquify(rev_it.locate()), (std::vector<uint64_t>{8}));
        EXPECT_RANGE_EQ(rev_it.path_label(this->rev_text1), this->pattern4);    // "GATGG"
        EXPECT_FALSE(rev_it.cycle_back());
    }

    {
        auto it = bi_fm.cursor();
        EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text, 3, 7))); // "GTAG"
        EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{3, 10}));

        auto rev_it = it.to_rev_cursor(); // text "CGATGCAGGATGGCA"
    #ifndef NDEBUG
        EXPECT_DEATH(rev_it.cycle_back(), "");
    #endif
        EXPECT_TRUE(rev_it.extend_right()); // "CGTAG" resp. "GATGC"
        EXPECT_EQ(seqan3::uniquify(rev_it.locate()), (std::vector<uint64_t>{1}));
        EXPECT_RANGE_EQ(rev_it.path_label(this->rev_text1), this->pattern3);    // "GATGC"
        EXPECT_TRUE(rev_it.cycle_back()); // "GGTAG" resp. "GATGG"
        EXPECT_EQ(seqan3::uniquify(rev_it.locate()), (std::vector<uint64_t>{8}));
        EXPECT_RANGE_EQ(rev_it.path_label(this->rev_text1), this->pattern4);    // "GATGG"
    }
}

REGISTER_TYPED_TEST_SUITE_P(bi_fm_index_cursor_test, cursor, extend, extend_char, extend_range, extend_and_cycle,
                            extend_range_and_cycle, to_fwd_cursor, to_rev_cursor);
