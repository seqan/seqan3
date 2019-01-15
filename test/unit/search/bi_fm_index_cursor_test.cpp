// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include "helper.hpp"

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/bi_fm_index_cursor.hpp>

#include <gtest/gtest.h>

using namespace seqan3;

template <typename T>
class bi_fm_index_cursor_test : public ::testing::Test
{};

using bi_fm_index_cursor_types = ::testing::Types<bi_fm_index_cursor<bi_fm_index<std::vector<dna4>>>>;

TYPED_TEST_CASE(bi_fm_index_cursor_test, bi_fm_index_cursor_types);

TYPED_TEST(bi_fm_index_cursor_test, begin)
{
    using text_t = typename TypeParam::index_type::text_type;
    text_t text{"AACGATCGGA"_dna4};
    auto rev_text = ranges::view::reverse(text);
    // using rev_text_t = decltype(ranges::view::reverse(text_t{}));

    using fm_fwd_t = fm_index<text_t>;
    using fm_rev_t = fm_index<decltype(rev_text)>;

    typename TypeParam::index_type bi_fm{text};
    fm_fwd_t fm_fwd{text};
    fm_rev_t fm_rev{rev_text};

    TypeParam bi_cur = bi_fm.begin();
    EXPECT_EQ(uniquify(bi_cur.locate()), uniquify(bi_fm.fwd_begin().locate()));
    EXPECT_EQ(uniquify(bi_cur.locate()), uniquify(bi_fm.rev_begin().locate()));
}

TYPED_TEST(bi_fm_index_cursor_test, extend)
{
    typename TypeParam::index_type::text_type text{"ACGGTAGGACG"_dna4};
    typename TypeParam::index_type bi_fm{text};

    auto cur = bi_fm.begin();
    EXPECT_TRUE(cur.extend_right()); // "A"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{0, 5, 8}));
    EXPECT_TRUE(cur.extend_left()); // "GA"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{7}));
    EXPECT_TRUE(cur.extend_right()); // "GAC"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{7}));
    EXPECT_TRUE(cur.extend_right()); // "GACG"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{7}));
    EXPECT_FALSE(cur.extend_right()); // "GACG"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{7}));
    EXPECT_TRUE(cur.extend_left()); // "GGACG"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{6}));
}

TYPED_TEST(bi_fm_index_cursor_test, extend_char)
{
    typename TypeParam::index_type::text_type text{"ACGGTAGGACG"_dna4};
    typename TypeParam::index_type bi_fm{text};

    auto cur = bi_fm.begin();
    EXPECT_TRUE(cur.extend_left('G'_dna4)); // "G"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{2, 3, 6, 7, 10}));
    EXPECT_TRUE(cur.extend_left('C'_dna4)); // "CG"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{1, 9}));
    EXPECT_FALSE(cur.extend_left('C'_dna4)); // "CG"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{1, 9}));
    EXPECT_FALSE(cur.extend_left('G'_dna4)); // "CG"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{1, 9}));
    EXPECT_FALSE(cur.extend_right('T'_dna4)); // "CG"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{1, 9}));
    EXPECT_TRUE(cur.extend_right('G'_dna4)); // "CGG"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{1}));
    EXPECT_TRUE(cur.extend_right('T'_dna4)); // "CGGT"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{1}));
    EXPECT_TRUE(cur.extend_right('A'_dna4)); // "CGGTA"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{1}));
    EXPECT_TRUE(cur.extend_left('A'_dna4)); // "ACGGTA"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{0}));
    EXPECT_FALSE(cur.extend_left('A'_dna4)); // "ACGGTA"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{0}));
}

TYPED_TEST(bi_fm_index_cursor_test, extend_range)
{
    typename TypeParam::index_type::text_type text{"ACGGTAGGACG"_dna4};
    typename TypeParam::index_type bi_fm{text};

    auto cur = bi_fm.begin();
    EXPECT_FALSE(cur.extend_left("CAG"_dna4)); // ""
    // sentinel poscurion included
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
    EXPECT_TRUE(cur.extend_left("CG"_dna4)); // "CG"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{1, 9}));
    EXPECT_TRUE(cur.extend_right("GTA"_dna4)); // "CGGTA"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{1}));
    EXPECT_FALSE(cur.extend_left("TA"_dna4)); // "CGGTA"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{1}));
    EXPECT_TRUE(cur.extend_left("A"_dna4)); // "ACGGTA"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{0}));
}

TYPED_TEST(bi_fm_index_cursor_test, extend_and_cycle)
{
    typename TypeParam::index_type::text_type text{"ACGGTAGGACG"_dna4};
    typename TypeParam::index_type bi_fm{text};

    auto cur = bi_fm.begin();
    EXPECT_TRUE(cur.extend_right()); // "A"
#ifndef NDEBUG
    EXPECT_DEATH(cur.cycle_front(), "");
#endif
    EXPECT_TRUE(cur.extend_left()); // "GA"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{7}));
#ifndef NDEBUG
    EXPECT_DEATH(cur.cycle_back(), "");
#endif
    EXPECT_TRUE(cur.cycle_front()); // "TA"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{4}));
    EXPECT_FALSE(cur.cycle_front()); // "TA"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{4}));
}

TYPED_TEST(bi_fm_index_cursor_test, extend_range_and_cycle)
{
    typename TypeParam::index_type::text_type text{"ACGGTAGGACGTAG"_dna4};
    typename TypeParam::index_type bi_fm{text};

    auto cur = bi_fm.begin();
    EXPECT_TRUE(cur.extend_right("AC"_dna4)); // "AC"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{0, 8}));
#ifndef NDEBUG
    EXPECT_DEATH(cur.cycle_front(), "");
#endif
    EXPECT_TRUE(cur.cycle_back()); // "AG"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{5, 12}));
#ifndef NDEBUG
    EXPECT_DEATH(cur.cycle_front(), "");
#endif
    EXPECT_FALSE(cur.extend_left("TT"_dna4)); // "AG"
    EXPECT_TRUE(cur.extend_left("CGT"_dna4)); // "CGTAG"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{9}));
#ifndef NDEBUG
    EXPECT_DEATH(cur.cycle_back(), "");
#endif
    EXPECT_TRUE(cur.cycle_front()); // "GGTAG"
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{2}));
}

TYPED_TEST(bi_fm_index_cursor_test, to_fwd_cursor)
{
    typename TypeParam::index_type::text_type text{"ACGGTAGGACGTAGC"_dna4};
    typename TypeParam::index_type bi_fm{text};

    {
        auto cur = bi_fm.begin();
        EXPECT_TRUE(cur.extend_right("GTAGC"_dna4)); // "GTAGC"
        EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{10}));

        auto fwd_cur = cur.to_fwd_cursor();
        EXPECT_TRUE(fwd_cur.cycle_back()); // "GTAGG"
        EXPECT_EQ(uniquify(fwd_cur.locate()), (std::vector<uint64_t>{3}));
        EXPECT_TRUE(std::ranges::equal(*fwd_cur, "GTAGG"_dna4));
        EXPECT_FALSE(fwd_cur.cycle_back());
    }

    {
        auto cur = bi_fm.begin();
        EXPECT_TRUE(cur.extend_left("GTAG"_dna4)); // "GTAG"
        EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{3, 10}));

        auto fwd_cur = cur.to_fwd_cursor();
    #ifndef NDEBUG
        EXPECT_DEATH(fwd_cur.cycle_back(), "");
    #endif
        EXPECT_TRUE(fwd_cur.extend_right());
        EXPECT_EQ(uniquify(fwd_cur.locate()), (std::vector<uint64_t>{10}));
        EXPECT_TRUE(std::ranges::equal(*fwd_cur, "GTAGC"_dna4));
        EXPECT_TRUE(fwd_cur.cycle_back());
        EXPECT_EQ(uniquify(fwd_cur.locate()), (std::vector<uint64_t>{3}));
        EXPECT_TRUE(std::ranges::equal(*fwd_cur, "GTAGG"_dna4));
    }
}

TYPED_TEST(bi_fm_index_cursor_test, to_rev_cursor)
{
    typename TypeParam::index_type::text_type text{"ACGGTAGGACGTAGC"_dna4};
    typename TypeParam::index_type bi_fm{text};

    {
        auto cur = bi_fm.begin();
        EXPECT_TRUE(cur.extend_left("CGTAG"_dna4)); // "CGTAG"
        EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{9}));

        auto rev_cur = cur.to_rev_cursor(); // text "CGATGCAGGATGGCA"
        EXPECT_EQ(uniquify(rev_cur.locate()), (std::vector<uint64_t>{1}));
        EXPECT_TRUE(std::ranges::equal(*rev_cur, "GATGC"_dna4));
        EXPECT_TRUE(rev_cur.cycle_back()); // "GATGG"
        EXPECT_EQ(uniquify(rev_cur.locate()), (std::vector<uint64_t>{8}));
        EXPECT_TRUE(std::ranges::equal(*rev_cur, "GATGG"_dna4));
        EXPECT_FALSE(rev_cur.cycle_back());
    }

    {
        auto cur = bi_fm.begin();
        EXPECT_TRUE(cur.extend_right("GTAG"_dna4)); // "GTAG"
        EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{3, 10}));

        auto rev_cur = cur.to_rev_cursor(); // text "CGATGCAGGATGGCA"
    #ifndef NDEBUG
        EXPECT_DEATH(rev_cur.cycle_back(), "");
    #endif
        EXPECT_TRUE(rev_cur.extend_right()); // "CGTAG" resp. "GATGC"
        EXPECT_EQ(uniquify(rev_cur.locate()), (std::vector<uint64_t>{1}));
        EXPECT_TRUE(std::ranges::equal(*rev_cur, "GATGC"_dna4));
        EXPECT_TRUE(rev_cur.cycle_back()); // "GGTAG" resp. "GATGG"
        EXPECT_EQ(uniquify(rev_cur.locate()), (std::vector<uint64_t>{8}));
        EXPECT_TRUE(std::ranges::equal(*rev_cur, "GATGG"_dna4));
    }
}
