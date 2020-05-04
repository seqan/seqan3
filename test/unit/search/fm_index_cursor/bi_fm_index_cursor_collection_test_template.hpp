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
struct bi_fm_index_cursor_collection_test;

TYPED_TEST_SUITE_P(bi_fm_index_cursor_collection_test);

TYPED_TEST_P(bi_fm_index_cursor_collection_test, cursor)
{
    typename TypeParam::index_type bi_fm{this->text_col1};  // {"AACGATCGGA", "AACGATCGGA"}

    seqan3::bi_fm_index fm_fwd{this->text_col1};
    seqan3::bi_fm_index fm_rev{this->rev_text1};

    TypeParam bi_it = bi_fm.cursor();
    EXPECT_EQ(seqan3::uniquify(bi_it.locate()), seqan3::uniquify(bi_fm.fwd_begin().locate()));
    EXPECT_EQ(seqan3::uniquify(bi_it.locate()), seqan3::uniquify(bi_fm.rev_begin().locate()));
}

TYPED_TEST_P(bi_fm_index_cursor_collection_test, extend)
{
    typename TypeParam::index_type bi_fm{this->text_col2};  // {"ACGGTAGGACG", "TGCTACGATCC"}

    auto it = bi_fm.cursor();
    EXPECT_TRUE(it.extend_right()); // "A"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 0}, {0, 5}, {0, 8},
                                                                                         {1, 4}, {1, 7}}));
    EXPECT_TRUE(it.extend_left()); // "GA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 7}, {1,6}}));
    EXPECT_TRUE(it.extend_right()); // "GAC"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 7}}));
    EXPECT_TRUE(it.extend_right()); // "GACG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 7}}));
    EXPECT_FALSE(it.extend_right()); // "GACG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 7}}));
    EXPECT_TRUE(it.extend_left()); // "GGACG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 6}}));
}

TYPED_TEST_P(bi_fm_index_cursor_collection_test, extend_char)
{
    typename TypeParam::index_type bi_fm{this->text_col2};  // {"ACGGTAGGACG", "TGCTACGATCC"}

    auto it = bi_fm.cursor();
    EXPECT_TRUE(it.extend_left(this->text[2])); // "G"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 2}, {0, 3}, {0, 6},
                                                                                         {0, 7}, {0, 10},
                                                                                         {1, 1}, {1, 6}}));
    EXPECT_TRUE(it.extend_left(this->text[1])); // "CG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}, {0, 9}, {1, 5}}));
    EXPECT_FALSE(it.extend_left(this->text[1])); // "CG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}, {0, 9}, {1, 5}}));
    EXPECT_FALSE(it.extend_left(this->text[2])); // "CG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}, {0, 9}, {1, 5}}));
    EXPECT_FALSE(it.extend_right(this->text[4])); // "CG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}, {0, 9}, {1, 5}}));
    EXPECT_TRUE(it.extend_right(this->text[2])); // "CGG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}}));
    EXPECT_TRUE(it.extend_right(this->text[4])); // "CGGT"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}}));
    EXPECT_TRUE(it.extend_right(this->text[0])); // "CGGTA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}}));
    EXPECT_TRUE(it.extend_left(this->text[0])); // "ACGGTA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 0}}));
    EXPECT_FALSE(it.extend_left(this->text[0])); // "ACGGTA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 0}}));
}

TYPED_TEST_P(bi_fm_index_cursor_collection_test, extend_range)
{
    typename TypeParam::index_type bi_fm{this->text_col2};  // {"ACGGTAGGACG", "TGCTACGATCC"}

    auto it = bi_fm.cursor();
    EXPECT_FALSE(it.extend_left(this->pattern1)); // "" ("CAG")
    // sentinel and delimiter position included
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 0}, {0, 1}, {0, 2}, {0, 3},
                                                                                         {0, 4}, {0, 5}, {0, 6}, {0, 7},
                                                                                         {0, 8}, {0, 9}, {0, 10},
                                                                                         {0, 11},
                                                                                         {1, 0}, {1, 1}, {1, 2}, {1, 3},
                                                                                         {1, 4}, {1, 5}, {1, 6}, {1, 7},
                                                                                         {1, 8}, {1, 9}, {1, 10},
                                                                                         {1, 11}}));
    EXPECT_TRUE(it.extend_left(seqan3::views::slice(this->text, 1, 3))); // "CG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}, {0, 9}, {1, 5}}));
    EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text, 3, 6))); // "CGGTA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}}));
    EXPECT_FALSE(it.extend_left(seqan3::views::slice(this->text, 2, 6))); // "CGGTA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}}));
    EXPECT_TRUE(it.extend_left(seqan3::views::slice(this->text, 0, 1))); // "ACGGTA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 0}}));
}

TYPED_TEST_P(bi_fm_index_cursor_collection_test, extend_and_cycle)
{
    typename TypeParam::index_type bi_fm{this->text_col2};  // {"ACGGTAGGACG", "TGCTACGATCC"}

    auto it = bi_fm.cursor();
    EXPECT_TRUE(it.extend_right()); // "A"
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_front(), "");
#endif
    EXPECT_TRUE(it.extend_left()); // "GA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 7}, {1, 6}}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_back(), "");
#endif
    EXPECT_TRUE(it.cycle_front()); // "TA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 4}, {1, 3}}));
    EXPECT_FALSE(it.cycle_front()); // "TA"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 4}, {1, 3}}));
}

TYPED_TEST_P(bi_fm_index_cursor_collection_test, extend_range_and_cycle)
{
    typename TypeParam::index_type bi_fm{this->text_col3};  // {"ACGGTAGGACGTAG", "TGCTACGATCC"}

    auto it = bi_fm.cursor();
    EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text, 0, 2))); // "AC"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 0}, {0, 8}, {1, 4}}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_front(), "");
#endif
    EXPECT_TRUE(it.cycle_back()); // "AG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 5}, {0, 12}}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_front(), "");
#endif
    EXPECT_FALSE(it.extend_left(this->pattern2)); // "AG"
    EXPECT_TRUE(it.extend_left(seqan3::views::slice(this->text, 9, 12))); // "CGTAG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 9}}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_back(), "");
#endif
    EXPECT_TRUE(it.cycle_front()); // "GGTAG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 2}}));
}

TYPED_TEST_P(bi_fm_index_cursor_collection_test, to_fwd_cursor)
{
    typename TypeParam::index_type bi_fm{this->text_col4};  // {"ACGGTAGGACGTAGC", "TGCTACGATCC"}

    {
        auto it = bi_fm.cursor();
        EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text, 10, 15))); // "GTAGC"
        EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 10}}));

        auto fwd_it = it.to_fwd_cursor();
        EXPECT_TRUE(fwd_it.cycle_back()); // "GTAGG"
        EXPECT_EQ(seqan3::uniquify(fwd_it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 3}}));
        EXPECT_RANGE_EQ(fwd_it.path_label(this->text_col4), seqan3::views::slice(this->text, 3, 8));
        EXPECT_FALSE(fwd_it.cycle_back());
    }

    {
        auto it = bi_fm.cursor();
        EXPECT_TRUE(it.extend_left(seqan3::views::slice(this->text, 3, 7))); // "GTAG"
        EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 3}, {0, 10}}));

        auto fwd_it = it.to_fwd_cursor();
    #ifndef NDEBUG
        EXPECT_DEATH(fwd_it.cycle_back(), "");
    #endif
        EXPECT_TRUE(fwd_it.extend_right());
        EXPECT_EQ(seqan3::uniquify(fwd_it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 10}}));
        EXPECT_RANGE_EQ(fwd_it.path_label(this->text_col4), seqan3::views::slice(this->text, 10, 15)); // "GTAGC"
        EXPECT_TRUE(fwd_it.cycle_back());
        EXPECT_EQ(seqan3::uniquify(fwd_it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 3}}));
        EXPECT_RANGE_EQ(fwd_it.path_label(this->text_col4), seqan3::views::slice(this->text, 3, 8)); // "GTAGG"
    }
}

TYPED_TEST_P(bi_fm_index_cursor_collection_test, to_rev_cursor)
{
    typename TypeParam::index_type bi_fm{this->text_col4};  // {"ACGGTAGGACGTAGC", "TGCTACGATCC"}

    {
        auto it = bi_fm.cursor();
        EXPECT_TRUE(it.extend_left(seqan3::views::slice(this->text, 9, 14))); // "CGTAG"
        EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 9}}));

        auto rev_it = it.to_rev_cursor(); // text_col4 "CCTAGCATCGT|CGATGCAGGATGGCA"
        EXPECT_EQ(seqan3::uniquify(rev_it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{1, 1}}));
        EXPECT_RANGE_EQ(rev_it.path_label(this->rev_text2), this->pattern3);   //ATGCA
        EXPECT_TRUE(rev_it.cycle_back()); // "GATGG"
        EXPECT_EQ(seqan3::uniquify(rev_it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{1, 8}}));
        EXPECT_RANGE_EQ(rev_it.path_label(this->rev_text2), this->pattern4);  // "GATGG"
        EXPECT_FALSE(rev_it.cycle_back());
    }

    {
        auto it = bi_fm.cursor();
        EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text, 3, 7))); // "GTAG"
        EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 3}, {0, 10}}));

        auto rev_it = it.to_rev_cursor(); // text_col4 "CCTAGCATCGT|CGATGCAGGATGGCA"
    #ifndef NDEBUG
        EXPECT_DEATH(rev_it.cycle_back(), "");
    #endif
        EXPECT_TRUE(rev_it.extend_right()); // "CGTAG" resp. "GATGC"
        EXPECT_EQ(seqan3::uniquify(rev_it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{1, 1}}));
        EXPECT_RANGE_EQ(rev_it.path_label(this->rev_text2), this->pattern3);  // "GATGC"
        EXPECT_TRUE(rev_it.cycle_back()); // "GGTAG" resp. "GATGG"
        EXPECT_EQ(seqan3::uniquify(rev_it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{1, 8}}));
        EXPECT_RANGE_EQ(rev_it.path_label(this->rev_text2), this->pattern4);  // "GATGG"
    }
}

TYPED_TEST_P(bi_fm_index_cursor_collection_test, extend_const_char_pointer)
{
    using alphabet_type = typename TestFixture::alphabet_type;

    // test case for https://github.com/seqan/seqan3/issues/1473
    if constexpr (std::is_same<alphabet_type, char>::value)
    {
        typename TypeParam::index_type fm{this->text_col1}; // {"AACGATCGGA", "AACGATCGGA"}
        char const * cg = "CG";

        // extend_right()
        {
            TypeParam it1 = TypeParam(fm);
            TypeParam it2 = TypeParam(fm);

            it1.extend_right(cg);
            it2.extend_right(seqan3::views::slice(this->text1, 1, 3));      // "CG"

            EXPECT_RANGE_EQ(it1.locate(), it2.locate());
        }
        // extend_left()
        {
            auto it1 = TypeParam(fm);
            auto it2 = TypeParam(fm);

            it1.extend_left(cg);
            it2.extend_right(seqan3::views::slice(this->text1, 1, 3));      // "CG"

            EXPECT_RANGE_EQ(it1.locate(), it2.locate());
        }
    }
}

REGISTER_TYPED_TEST_SUITE_P(bi_fm_index_cursor_collection_test, cursor, extend, extend_char, extend_range,
                            extend_and_cycle, extend_range_and_cycle, to_fwd_cursor, to_rev_cursor, extend_const_char_pointer);
