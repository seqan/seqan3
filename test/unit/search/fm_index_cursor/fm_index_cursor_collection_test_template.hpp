/// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include "../helper.hpp"

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/std/algorithm>

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
class fm_index_cursor_collection_test : public ::testing::Test
{};

TYPED_TEST_CASE_P(fm_index_cursor_collection_test);

TYPED_TEST_P(fm_index_cursor_collection_test, ctr)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4, "ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // custom constructor
    TypeParam it0{fm};
    EXPECT_EQ(it0.query_length(), 0u);
    EXPECT_EQ(it0.locate().size(), fm.size());

    // default construction (does not initialize the cursor)
    [[maybe_unused]] TypeParam it1;

    // copy construction
    TypeParam it2{it0};
    EXPECT_EQ(it0, it2);

    // copy assignment
    TypeParam it3 = it0;
    EXPECT_EQ(it0, it3);

    // move construction
    TypeParam it4{std::move(it0)};
    EXPECT_EQ(it0, it4);

    // move assigment
    TypeParam it5 = std::move(it0);
    EXPECT_EQ(it0, it5);
}

TYPED_TEST_P(fm_index_cursor_collection_test, begin)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4, "ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // begin
    TypeParam it(fm);
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,0}, {0,1}, {0,2}, {0,3}, {0,4},
                                                                                 {0,5}, {0,6}, {1,0},{1,1}, {1,2},
                                                                                 {1,3}, {1,4}, {1,5}, {1,6}}));
                                                                                 // one sentinel position included
    EXPECT_EQ(it.query_length(), 0u);
    EXPECT_EQ(it.count(), 14u);
}

TYPED_TEST_P(fm_index_cursor_collection_test, extend_right_range)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4, "TGCGATCGA"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right(range)
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right("CG"_dna4));
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,1}, {0,4}, {1,2}, {1,6}}));
    EXPECT_EQ(it.query_length(), 2u);
    EXPECT_EQ(it.count(), 4u);

    EXPECT_TRUE(it.extend_right("A"_dna4));
    EXPECT_EQ(it.locate(), (std::vector<std::pair<uint64_t, uint64_t>>{{0,1}, {1,2}, {1,6}}));
    EXPECT_EQ(it.query_length(), 3u);
    EXPECT_EQ(it.count(), 3u);

    // unsuccessful extend_right(range), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.extend_right("A"_dna4));
    EXPECT_EQ(it, it_cpy);

    // extend_right(empty range)
    it_cpy = it;
    EXPECT_TRUE(it.extend_right(""_dna4));
    EXPECT_EQ(it, it_cpy);
}

TYPED_TEST_P(fm_index_cursor_collection_test, extend_right_range_empty_text)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4, ""_dna4, ""_dna4, "TGCGATCGA"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right(range)
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right("CG"_dna4));
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,1}, {0,4}, {3,2}, {3,6}}));
    EXPECT_EQ(it.query_length(), 2u);
    EXPECT_EQ(it.count(), 4u);

    EXPECT_TRUE(it.extend_right("A"_dna4));
    EXPECT_EQ(it.locate(), (std::vector<std::pair<uint64_t, uint64_t>>{{0,1}, {3,2}, {3,6}}));
    EXPECT_EQ(it.query_length(), 3u);
    EXPECT_EQ(it.count(), 3u);

    // unsuccessful extend_right(range), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.extend_right("A"_dna4));
    EXPECT_EQ(it, it_cpy);

    // extend_right(empty range)
    it_cpy = it;
    EXPECT_TRUE(it.extend_right(""_dna4));
    EXPECT_EQ(it, it_cpy);
}

// TODO: doesn't work with the current structure of typed tests
// TYPED_TEST_P(fm_index_cursor_collection_test, extend_right_convertible_range)
// {
//     typename TypeParam::index_type::text_type text{"ANGACGNN"_dna5};
//     typename TypeParam::index_type fm{text};
//
//     // successful extend_right(range) using a different alphabet
//     TypeParam it(fm);
//     EXPECT_TRUE(it.extend_right("GA"_dna4));
//     EXPECT_EQ(it.locate(), (std::vector<uint64_t>{2}));
//     EXPECT_EQ(it.query_length(), 2);
// }

TYPED_TEST_P(fm_index_cursor_collection_test, extend_right_char)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4, "TGCGATCGA"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right(char)
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right("A"_dna4));
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,0}, {0,3}, {1,4}, {1,8}}));
    EXPECT_EQ(it.query_length(), 1u);

    EXPECT_TRUE(it.extend_right("C"_dna4));
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,0}, {0,3}}));
    EXPECT_EQ(it.query_length(), 2u);

    // unsuccessful extend_right(char), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.extend_right("C"_dna4));
    EXPECT_EQ(it, it_cpy);
}

// TODO: doesn't work with the current structure of typed tests
// TYPED_TEST_P(fm_index_cursor_collection_test, extend_right_convertible_char)
// {
//     typename TypeParam::index_type::text_type text{"ANGACGNN"_dna5};
//     typename TypeParam::index_type fm{text};
//
//     // successful extend_right(char) using a different alphabet
//     TypeParam it(fm);
//     EXPECT_TRUE(it.extend_right("A"_dna4));
//     EXPECT_EQ(uniquify(it.locate()), (std::vector<uint64_t>{0, 3}));
//     EXPECT_EQ(it.query_length(), 1);
// }

TYPED_TEST_P(fm_index_cursor_collection_test, extend_right_range_and_cycle)
{
    typename TypeParam::index_type::text_type text{"ACGAACGC"_dna4, "TACGATCGA"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right() and cycle_back()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right("ACGA"_dna4));
    EXPECT_EQ(it.locate(), (std::vector<std::pair<uint64_t, uint64_t>>{{0,0}, {1,1}}));
    EXPECT_EQ(it.query_length(), 4u);

    EXPECT_TRUE(it.cycle_back());
    EXPECT_EQ(it.locate(), (std::vector<std::pair<uint64_t, uint64_t>>{{0,4}}));
    EXPECT_EQ(it.query_length(), 4u);
}

TYPED_TEST_P(fm_index_cursor_collection_test, extend_right_char_and_cycle)
{
    typename TypeParam::index_type::text_type text{"ACGAACGC"_dna4, "TGCGATCGA"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right() and cycle_back()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right("A"_dna4));
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,0}, {0,3}, {0,4}, {1,4}, {1,8}}));
    EXPECT_EQ(it.query_length(), 1u);

    EXPECT_TRUE(it.cycle_back());
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,1}, {0,5}, {0,7}, {1,2}, {1,6}}));
    EXPECT_EQ(it.query_length(), 1u);
}

TYPED_TEST_P(fm_index_cursor_collection_test, extend_right_and_cycle)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4, "TGCGATCGA"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right() and cycle_back()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right());
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,0}, {0,3}, {1,4}, {1,8}}));
    EXPECT_EQ(it.query_length(), 1u);

    EXPECT_TRUE(it.cycle_back());
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,1}, {0,4}, {1,2}, {1,6}}));
    EXPECT_EQ(it.query_length(), 1u);

    EXPECT_TRUE(it.extend_right());
    EXPECT_EQ(uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,1}, {0,4}, {1,2}, {1,6}}));
    EXPECT_EQ(it.query_length(), 2u);

    // unsuccessful cycle_back(), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.cycle_back());
    EXPECT_EQ(it, it_cpy);

    // unsuccessful extend_right(), it remains untouched
    it = TypeParam(fm);
    EXPECT_TRUE(it.extend_right("ACGACG"_dna4));
    it_cpy = it;
    EXPECT_FALSE(it.extend_right());
    EXPECT_EQ(it, it_cpy);

    // cycle_back() cannot be called on the root node
    it = TypeParam(fm);
#ifndef NDEBUG
    EXPECT_DEATH(it.cycle_back(), "");
#endif
    EXPECT_EQ(it, TypeParam(fm));
}

TYPED_TEST_P(fm_index_cursor_collection_test, query)
{
    typename TypeParam::index_type::text_type text{"ACGACG"_dna4, "TGCGATCGA"_dna4};
    typename TypeParam::index_type fm{text};

    // query()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right("ACG"_dna4));
    EXPECT_TRUE(std::ranges::equal(*it, "ACG"_dna4));
    EXPECT_TRUE(std::ranges::equal(it.query(), "ACG"_dna4));
}

// TODO: test last_char()
TYPED_TEST_P(fm_index_cursor_collection_test, incomplete_alphabet)
{
    // search a char that does not occur in the text (higher rank than largest char occurring in text)
    {
        typename TypeParam::index_type::text_type text{"ACGACG"_dna4, "ACGACG"_dna4};
        typename TypeParam::index_type fm{text};
        TypeParam it = TypeParam(fm);
        EXPECT_FALSE(it.extend_right("T"_dna4));
        EXPECT_EQ(it, TypeParam(fm));
    }

    // search a char that does not occur in the text (smaller rank than smallest char occurring in text)
    {
        typename TypeParam::index_type::text_type text{"CGTCGT"_dna4, "CGTCGT"_dna4};
        typename TypeParam::index_type fm{text};
        TypeParam it = TypeParam(fm);
        EXPECT_FALSE(it.extend_right("A"_dna4));
        EXPECT_EQ(it, TypeParam(fm));
    }

    // search a char that does not occur in the text
    // (some rank that is neither the smallest nor the highest occurring in text)
    {
        typename TypeParam::index_type::text_type text{"ATATAT"_dna4, "ATATAT"_dna4};
        typename TypeParam::index_type fm{text};
        TypeParam it = TypeParam(fm);
        EXPECT_FALSE(it.extend_right("C"_dna4));
        EXPECT_FALSE(it.extend_right("G"_dna4));
        EXPECT_FALSE(it.extend_right("ACGT"_dna4));
        EXPECT_FALSE(it.extend_right("G"_dna4));
        EXPECT_EQ(it, TypeParam(fm));

        EXPECT_TRUE(it.extend_right("A"_dna4));
        EXPECT_TRUE(it.cycle_back());
        EXPECT_TRUE(std::ranges::equal(it.query(), "T"_dna4));
    }
}

TYPED_TEST_P(fm_index_cursor_collection_test, lazy_locate)
{
    typename TypeParam::index_type::text_type text{"ACGTACGT"_dna4, "TGCGATACGA"_dna4};
    typename TypeParam::index_type fm{text};

    TypeParam it = TypeParam(fm);
    it.extend_right("ACG"_dna4);

    EXPECT_TRUE(std::ranges::equal(it.locate(), it.lazy_locate()));
}

TYPED_TEST_P(fm_index_cursor_collection_test, concept_check)
{
    EXPECT_TRUE(FmIndexCursor<TypeParam>);
}

REGISTER_TYPED_TEST_CASE_P(fm_index_cursor_collection_test, ctr, begin, extend_right_range,
                           extend_right_range_empty_text, extend_right_char, extend_right_range_and_cycle,
                           extend_right_char_and_cycle, extend_right_and_cycle, query, incomplete_alphabet, lazy_locate,
                           concept_check);
