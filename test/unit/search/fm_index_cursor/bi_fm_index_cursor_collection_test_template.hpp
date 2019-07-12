/// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include "../helper.hpp"

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/bi_fm_index_cursor.hpp>
#include <seqan3/std/algorithm>

#include <gtest/gtest.h>

using namespace seqan3;

template <typename T>
class bi_fm_index_cursor_collection_test : public ::testing::Test
{};

TYPED_TEST_CASE_P(bi_fm_index_cursor_collection_test);

TYPED_TEST_P(bi_fm_index_cursor_collection_test, begin)
{
    std::vector<std::vector<dna4>> text{"AACGATCGGA"_dna4, "AACGATCGGA"_dna4};
    auto rev_text = text | view::deep{std::view::reverse} | view::deep{view::persist};

    typename TypeParam::index_type bi_fm{text};
    bi_fm_index fm_fwd{text};
    bi_fm_index fm_rev{rev_text};

    TypeParam bi_it = bi_fm.begin();
    EXPECT_EQ(uniquify(bi_it.locate()), uniquify(bi_fm.fwd_begin().locate()));
    EXPECT_EQ(uniquify(bi_it.locate()), uniquify(bi_fm.rev_begin().locate()));
}

TYPED_TEST_P(bi_fm_index_cursor_collection_test, extend)
{
    std::vector<std::vector<dna4>> text{"ACGGTAGGACG"_dna4, "TGCTACGATCC"_dna4};
    typename TypeParam::index_type bi_fm{text};

    auto it = bi_fm.begin();
    EXPECT_TRUE(it.extend_right()); // "A"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 0}, {0, 5}, {0, 8},
                                                                                 {1, 4}, {1, 7}}));
    EXPECT_TRUE(it.extend_left()); // "GA"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 7}, {1,6}}));
    EXPECT_TRUE(it.extend_right()); // "GAC"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 7}}));
    EXPECT_TRUE(it.extend_right()); // "GACG"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 7}}));
    EXPECT_FALSE(it.extend_right()); // "GACG"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 7}}));
    EXPECT_TRUE(it.extend_left()); // "GGACG"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 6}}));
}

TYPED_TEST_P(bi_fm_index_cursor_collection_test, extend_char)
{
    std::vector<std::vector<dna4>> text{"ACGGTAGGACG"_dna4, "TGCTACGATCC"_dna4};
    typename TypeParam::index_type bi_fm{text};

    auto it = bi_fm.begin();
    EXPECT_TRUE(it.extend_left('G'_dna4)); // "G"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 2}, {0, 3}, {0, 6},
                                                                                 {0, 7}, {0, 10},
                                                                                 {1, 1}, {1, 6}}));
    EXPECT_TRUE(it.extend_left('C'_dna4)); // "CG"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}, {0, 9}, {1, 5}}));
    EXPECT_FALSE(it.extend_left('C'_dna4)); // "CG"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}, {0, 9}, {1, 5}}));
    EXPECT_FALSE(it.extend_left('G'_dna4)); // "CG"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}, {0, 9}, {1, 5}}));
    EXPECT_FALSE(it.extend_right('T'_dna4)); // "CG"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}, {0, 9}, {1, 5}}));
    EXPECT_TRUE(it.extend_right('G'_dna4)); // "CGG"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}}));
    EXPECT_TRUE(it.extend_right('T'_dna4)); // "CGGT"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}}));
    EXPECT_TRUE(it.extend_right('A'_dna4)); // "CGGTA"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}}));
    EXPECT_TRUE(it.extend_left('A'_dna4)); // "ACGGTA"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 0}}));
    EXPECT_FALSE(it.extend_left('A'_dna4)); // "ACGGTA"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 0}}));
}

TYPED_TEST_P(bi_fm_index_cursor_collection_test, extend_range)
{
    std::vector<std::vector<dna4>> text{"ACGGTAGGACG"_dna4, "TGCTACGATCC"_dna4};
    typename TypeParam::index_type bi_fm{text};

    auto it = bi_fm.begin();
    EXPECT_FALSE(it.extend_left("CAG"_dna4)); // ""
    // sentinel and delimiter position included
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4},
                                                                                 {0, 5}, {0, 6}, {0, 7}, {0, 8}, {0, 9},
                                                                                 {0, 10}, {0, 11},
                                                                                 {1, 0}, {1, 1}, {1, 2}, {1, 3}, {1, 4},
                                                                                 {1, 5}, {1, 6}, {1, 7}, {1, 8}, {1, 9},
                                                                                 {1, 10}, {1, 11}}));
    EXPECT_TRUE(it.extend_left("CG"_dna4)); // "CG"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}, {0, 9}, {1, 5}}));
    EXPECT_TRUE(it.extend_right("GTA"_dna4)); // "CGGTA"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}}));
    EXPECT_FALSE(it.extend_left("TA"_dna4)); // "CGGTA"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 1}}));
    EXPECT_TRUE(it.extend_left("A"_dna4)); // "ACGGTA"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 0}}));
}

TYPED_TEST_P(bi_fm_index_cursor_collection_test, extend_and_cycle)
{
    std::vector<std::vector<dna4>> text{"ACGGTAGGACG"_dna4, "TGCTACGATCC"_dna4};
    typename TypeParam::index_type bi_fm{text};

    auto it = bi_fm.begin();
    EXPECT_TRUE(it.extend_right()); // "A"
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_front(), "");
#endif
    EXPECT_TRUE(it.extend_left()); // "GA"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 7}, {1, 6}}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_back(), "");
#endif
    EXPECT_TRUE(it.cycle_front()); // "TA"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 4}, {1, 3}}));
    EXPECT_FALSE(it.cycle_front()); // "TA"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 4}, {1, 3}}));
}

TYPED_TEST_P(bi_fm_index_cursor_collection_test, extend_range_and_cycle)
{
    std::vector<std::vector<dna4>> text{"ACGGTAGGACGTAG"_dna4, "TGCTACGATCC"_dna4};
    typename TypeParam::index_type bi_fm{text};

    auto it = bi_fm.begin();
    EXPECT_TRUE(it.extend_right("AC"_dna4)); // "AC"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 0}, {0, 8}, {1, 4}}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_front(), "");
#endif
    EXPECT_TRUE(it.cycle_back()); // "AG"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 5}, {0, 12}}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_front(), "");
#endif
    EXPECT_FALSE(it.extend_left("TT"_dna4)); // "AG"
    EXPECT_TRUE(it.extend_left("CGT"_dna4)); // "CGTAG"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 9}}));
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_back(), "");
#endif
    EXPECT_TRUE(it.cycle_front()); // "GGTAG"
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 2}}));
}

TYPED_TEST_P(bi_fm_index_cursor_collection_test, to_fwd_cursor)
{
    std::vector<std::vector<dna4>> text{"ACGGTAGGACGTAGC"_dna4, "TGCTACGATCC"_dna4};
    typename TypeParam::index_type bi_fm{text};

    {
        auto it = bi_fm.begin();
        EXPECT_TRUE(it.extend_right("GTAGC"_dna4)); // "GTAGC"
        EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 10}}));

        auto fwd_it = it.to_fwd_cursor();
        EXPECT_TRUE(fwd_it.cycle_back()); // "GTAGG"
        EXPECT_EQ(uniquify(fwd_it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 3}}));
        EXPECT_TRUE(std::ranges::equal(fwd_it.path_label(text), "GTAGG"_dna4));
        EXPECT_FALSE(fwd_it.cycle_back());
    }

    {
        auto it = bi_fm.begin();
        EXPECT_TRUE(it.extend_left("GTAG"_dna4)); // "GTAG"
        EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 3}, {0, 10}}));

        auto fwd_it = it.to_fwd_cursor();
    #ifndef NDEBUG
        EXPECT_DEATH(fwd_it.cycle_back(), "");
    #endif
        EXPECT_TRUE(fwd_it.extend_right());
        EXPECT_EQ(uniquify(fwd_it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 10}}));
        EXPECT_TRUE(std::ranges::equal(fwd_it.path_label(text), "GTAGC"_dna4));
        EXPECT_TRUE(fwd_it.cycle_back());
        EXPECT_EQ(uniquify(fwd_it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 3}}));
        EXPECT_TRUE(std::ranges::equal(fwd_it.path_label(text), "GTAGG"_dna4));
    }
}

TYPED_TEST_P(bi_fm_index_cursor_collection_test, to_rev_cursor)
{
    std::vector<std::vector<dna4>> text{"ACGGTAGGACGTAGC"_dna4, "TGCTACGATCC"_dna4};
    std::vector<std::vector<dna4>> rev_text{text | view::deep{std::view::reverse} | std::view::reverse};
    typename TypeParam::index_type bi_fm{text};

    {
        auto it = bi_fm.begin();
        EXPECT_TRUE(it.extend_left("CGTAG"_dna4)); // "CGTAG"
        EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 9}}));

        auto rev_it = it.to_rev_cursor(); // text "CCTAGCATCGT|CGATGCAGGATGGCA"
        EXPECT_EQ(uniquify(rev_it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{1, 1}}));
        EXPECT_TRUE(std::ranges::equal(rev_it.path_label(rev_text), "GATGC"_dna4));   //ATGCA
        EXPECT_TRUE(rev_it.cycle_back()); // "GATGG"
        EXPECT_EQ(uniquify(rev_it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{1, 8}}));
        EXPECT_TRUE(std::ranges::equal(rev_it.path_label(rev_text), "GATGG"_dna4));
        EXPECT_FALSE(rev_it.cycle_back());
    }

    {
        auto it = bi_fm.begin();
        EXPECT_TRUE(it.extend_right("GTAG"_dna4)); // "GTAG"
        EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0, 3}, {0, 10}}));

        auto rev_it = it.to_rev_cursor(); // text "CCTAGCATCGT|CGATGCAGGATGGCA"
    #ifndef NDEBUG
        EXPECT_DEATH(rev_it.cycle_back(), "");
    #endif
        EXPECT_TRUE(rev_it.extend_right()); // "CGTAG" resp. "GATGC"
        EXPECT_EQ(uniquify(rev_it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{1, 1}}));
        EXPECT_TRUE(std::ranges::equal(rev_it.path_label(rev_text), "GATGC"_dna4));
        EXPECT_TRUE(rev_it.cycle_back()); // "GGTAG" resp. "GATGG"
        EXPECT_EQ(uniquify(rev_it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{1, 8}}));
        EXPECT_TRUE(std::ranges::equal(rev_it.path_label(rev_text), "GATGG"_dna4));
    }
}

REGISTER_TYPED_TEST_CASE_P(bi_fm_index_cursor_collection_test, begin, extend, extend_char, extend_range,
                           extend_and_cycle, extend_range_and_cycle, to_fwd_cursor, to_rev_cursor);
