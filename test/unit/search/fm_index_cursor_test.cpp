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
#include <seqan3/search/fm_index/all.hpp>


#include <gtest/gtest.h>

using namespace seqan3;

struct fm_index_byte_alphabet_traits
{
    using sdsl_index_type = sdsl::csa_wt<
        sdsl::wt_blcd<
            sdsl::bit_vector,
            sdsl::rank_support_v<>,
            sdsl::select_support_scan<>,
            sdsl::select_support_scan<0>
        >,
        16,
        10000000,
        sdsl::sa_order_sa_sampling<>,
        sdsl::isa_sampling<>,
        sdsl::byte_alphabet
    >;
};

struct bi_fm_index_byte_alphabet_traits
{
   using fm_index_traits = fm_index_byte_alphabet_traits;
   using rev_fm_index_traits = fm_index_byte_alphabet_traits;
};

template <typename T>
class fm_index_cursor_test : public ::testing::Test
{};

using fm_index_cursor_types = ::testing::Types<
        fm_index_cursor<fm_index<std::vector<dna4>, fm_index_default_traits>>,
        fm_index_cursor<fm_index<std::vector<dna4>, fm_index_byte_alphabet_traits>>,
        bi_fm_index_cursor<bi_fm_index<std::vector<dna4>, bi_fm_index_default_traits>>,
        bi_fm_index_cursor<bi_fm_index<std::vector<dna4>, bi_fm_index_byte_alphabet_traits>>>;

TYPED_TEST_CASE(fm_index_cursor_test, fm_index_cursor_types);

TYPED_TEST(fm_index_cursor_test, ctr)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // custom constructor
    TypeParam cur0{fm};
    EXPECT_EQ(cur0.query_length(), 0u);
    EXPECT_EQ(cur0.locate().size(), fm.size());

    // default construction (does not initialize the cursor)
    TypeParam cur1;

    // copy construction
    TypeParam cur2{cur0};
    EXPECT_EQ(cur0, cur2);

    // copy assignment
    TypeParam cur3 = cur0;
    EXPECT_EQ(cur0, cur3);

    // move construction
    TypeParam cur4{std::move(cur0)};
    EXPECT_EQ(cur0, cur4);

    // move assigment
    TypeParam cur5 = std::move(cur0);
    EXPECT_EQ(cur0, cur5);
}

TYPED_TEST(fm_index_cursor_test, begin)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // begin
    TypeParam cur(fm);
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{0, 1, 2, 3, 4, 5, 6})); // sentinel position included
    EXPECT_EQ(cur.query_length(), 0u);
    EXPECT_EQ(cur.count(), 7u);
}

TYPED_TEST(fm_index_cursor_test, extend_right_range)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right(range)
    TypeParam cur(fm);
    EXPECT_TRUE(cur.extend_right("CG"_dna4));
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{1, 4}));
    EXPECT_EQ(cur.query_length(), 2u);
    EXPECT_EQ(cur.count(), 2u);

    EXPECT_TRUE(cur.extend_right("A"_dna4));
    EXPECT_EQ(cur.locate(), (std::vector<uint64_t>{1}));
    EXPECT_EQ(cur.query_length(), 3u);
    EXPECT_EQ(cur.count(), 1u);

    // unsuccessful extend_right(range), cur remains untouched
    TypeParam cur_cpy = cur;
    EXPECT_FALSE(cur.extend_right("A"_dna4));
    EXPECT_EQ(cur, cur_cpy);

    // extend_right(empty range)
    cur_cpy = cur;
    EXPECT_TRUE(cur.extend_right(""_dna4));
    EXPECT_EQ(cur, cur_cpy);
}

// TODO: doesn't work with the current structure of typed tests
// TYPED_TEST(fm_index_cursor_test, extend_right_convertible_range)
// {
//     typename TypeParam::index_type::text_type text{"ANGACGNN"_dna5};
//     typename TypeParam::index_type fm{text};
//
//     // successful extend_right(range) using a different alphabet
//     TypeParam cur(fm);
//     EXPECT_TRUE(cur.extend_right("GA"_dna4));
//     EXPECT_EQ(cur.locate(), (std::vector<uint64_t>{2}));
//     EXPECT_EQ(cur.query_length(), 2u);
// }

TYPED_TEST(fm_index_cursor_test, extend_right_char)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right(char)
    TypeParam cur(fm);
    EXPECT_TRUE(cur.extend_right('A'_dna4));
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{0, 3}));
    EXPECT_EQ(cur.query_length(), 1u);

    EXPECT_TRUE(cur.extend_right('C'_dna4));
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{0, 3}));
    EXPECT_EQ(cur.query_length(), 2u);

    // unsuccessful extend_right(char), cur remains untouched
    TypeParam cur_cpy = cur;
    EXPECT_FALSE(cur.extend_right('C'_dna4));
    EXPECT_EQ(cur, cur_cpy);
}

// TODO: doesn't work with the current structure of typed tests
// TYPED_TEST(fm_index_cursor_test, extend_right_convertible_char)
// {
//     typename TypeParam::index_type::text_type text{"ANGACGNN"_dna5};
//     typename TypeParam::index_type fm{text};
//
//     // successful extend_right(char) using a different alphabet
//     TypeParam cur(fm);
//     EXPECT_TRUE(cur.extend_right('A'_dna4));
//     EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{0, 3}));
//     EXPECT_EQ(cur.query_length(), 1u);
// }

TYPED_TEST(fm_index_cursor_test, extend_right_range_and_cycle)
{
    typename TypeParam::index_type::text_type text{"ACGAACGC"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right() and cycle_back()
    TypeParam cur(fm);
    EXPECT_TRUE(cur.extend_right("ACGA"_dna4));
    EXPECT_EQ(cur.locate(), (std::vector<uint64_t>{0}));
    EXPECT_EQ(cur.query_length(), 4u);

    EXPECT_TRUE(cur.cycle_back());
    EXPECT_EQ(cur.locate(), (std::vector<uint64_t>{4}));
    EXPECT_EQ(cur.query_length(), 4u);
}

TYPED_TEST(fm_index_cursor_test, extend_right_char_and_cycle)
{
    typename TypeParam::index_type::text_type text{"ACGAACGC"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right() and cycle_back()
    TypeParam cur(fm);
    EXPECT_TRUE(cur.extend_right('A'_dna4));
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{0, 3, 4}));
    EXPECT_EQ(cur.query_length(), 1u);

    EXPECT_TRUE(cur.cycle_back());
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{1, 5, 7}));
    EXPECT_EQ(cur.query_length(), 1u);
}

TYPED_TEST(fm_index_cursor_test, extend_right_and_cycle)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right() and cycle_back()
    TypeParam cur(fm);
    EXPECT_TRUE(cur.extend_right());
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{0, 3}));
    EXPECT_EQ(cur.query_length(), 1u);

    EXPECT_TRUE(cur.cycle_back());
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{1, 4}));
    EXPECT_EQ(cur.query_length(), 1u);

    EXPECT_TRUE(cur.extend_right());
    EXPECT_EQ(uniquify(cur.locate()), (std::vector<uint64_t>{1, 4}));
    EXPECT_EQ(cur.query_length(), 2u);

    // unsuccessful cycle_back(), cur remains untouched
    TypeParam cur_cpy = cur;
    EXPECT_FALSE(cur.cycle_back());
    EXPECT_EQ(cur, cur_cpy);

    // unsuccessful extend_right(), cur remains untouched
    cur = TypeParam(fm);
    EXPECT_TRUE(cur.extend_right("GACG"_dna4));
    cur_cpy = cur;
    EXPECT_FALSE(cur.extend_right());
    EXPECT_EQ(cur, cur_cpy);

    // cycle_back() cannot be called on the root node
    cur = TypeParam(fm);
#ifndef NDEBUG
    EXPECT_DEATH(cur.cycle_back(), "");
#endif
    EXPECT_EQ(cur, TypeParam(fm));
}

TYPED_TEST(fm_index_cursor_test, query)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // query()
    TypeParam cur(fm);
    EXPECT_TRUE(cur.extend_right("ACG"_dna4));
    EXPECT_TRUE(std::ranges::equal(*cur, "ACG"_dna4));
    EXPECT_TRUE(std::ranges::equal(cur.query(), "ACG"_dna4));
}

// TODO: test last_char()
TYPED_TEST(fm_index_cursor_test, incomplete_alphabet)
{
    // search a char that does not occur in the text (higher rank than largest char occurring in text)
    {
        typename TypeParam::index_type::text_type text{"ACGACG"_dna4};
        typename TypeParam::index_type fm{text};
        TypeParam cur = TypeParam(fm);
        EXPECT_FALSE(cur.extend_right('T'_dna4));
        EXPECT_EQ(cur, TypeParam(fm));
    }

    // search a char that does not occur in the text (smaller rank than smallest char occurring in text)
    {
        typename TypeParam::index_type::text_type text{"CGTCGT"_dna4};
        typename TypeParam::index_type fm{text};
        TypeParam cur = TypeParam(fm);
        EXPECT_FALSE(cur.extend_right('A'_dna4));
        EXPECT_EQ(cur, TypeParam(fm));
    }

    // search a char that does not occur in the text
    // (some rank that is neither the smallest nor the highest occurring in text)
    {
        typename TypeParam::index_type::text_type text{"ATATAT"_dna4};
        typename TypeParam::index_type fm{text};
        TypeParam cur = TypeParam(fm);
        EXPECT_FALSE(cur.extend_right('C'_dna4));
        EXPECT_FALSE(cur.extend_right('G'_dna4));
        EXPECT_FALSE(cur.extend_right("ACGT"_dna4));
        EXPECT_FALSE(cur.extend_right("G"_dna4));
        EXPECT_EQ(cur, TypeParam(fm));

        EXPECT_TRUE(cur.extend_right('A'_dna4));
        EXPECT_TRUE(cur.cycle_back());
        EXPECT_TRUE(std::ranges::equal(cur.query(), "T"_dna4));
    }
}

TYPED_TEST(fm_index_cursor_test, lazy_locate)
{
    typename TypeParam::index_type::text_type text{"ACGTACGT"_dna4};
    typename TypeParam::index_type fm{text};

    TypeParam cur = TypeParam(fm);
    cur.extend_right("ACG"_dna4);

    EXPECT_TRUE(std::ranges::equal(cur.locate(), cur.lazy_locate()));
}

TEST(fm_index, concepts)
{
    EXPECT_TRUE(FmIndexCursor<fm_index_cursor<fm_index<std::vector<dna4>>>>);
    EXPECT_TRUE(FmIndexCursor<bi_fm_index_cursor<bi_fm_index<std::vector<dna4>>>>);
    EXPECT_TRUE(BiFmIndexCursor<bi_fm_index_cursor<bi_fm_index<std::vector<dna4>>>>);
}
