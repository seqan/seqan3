/// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>
#include <type_traits>

#include <seqan3/range/views/slice.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/std/algorithm>

#include "../helper.hpp"

using sdsl_byte_index_type = sdsl::csa_wt<
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

template <typename T>
struct fm_index_cursor_collection_test;

TYPED_TEST_SUITE_P(fm_index_cursor_collection_test);

TYPED_TEST_P(fm_index_cursor_collection_test, ctr)
{
    typename TypeParam::index_type fm{this->text_col1}; // {"ACGACG", "ACGACG"}

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
    typename TypeParam::index_type fm{this->text_col1}; // {"ACGACG", "ACGACG"}

    // begin
    TypeParam it(fm);
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,0}, {0,1}, {0,2}, {0,3},
                                                                                         {0,4}, {0,5}, {0,6},
                                                                                         {1,0},{1,1}, {1,2}, {1,3},
                                                                                         {1,4}, {1,5}, {1,6}}));
                                                                                       // one sentinel position included
    EXPECT_EQ(it.query_length(), 0u);
    EXPECT_EQ(it.count(), 14u);
}

TYPED_TEST_P(fm_index_cursor_collection_test, extend_right_range)
{
    typename TypeParam::index_type fm{this->text_col2}; // {"ACGACG", "TGCGATCGA"}

    // successful extend_right(range)
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text1, 1, 3)));   // "CG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,1}, {0,4}, {1,2}, {1,6}}));
    EXPECT_EQ(it.query_length(), 2u);
    EXPECT_EQ(it.count(), 4u);

    EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text1, 0, 1)));   // "A"
    EXPECT_EQ(it.locate(), (std::vector<std::pair<uint64_t, uint64_t>>{{0,1}, {1,2}, {1,6}}));
    EXPECT_EQ(it.query_length(), 3u);
    EXPECT_EQ(it.count(), 3u);

    // unsuccessful extend_right(range), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.extend_right(seqan3::views::slice(this->text1, 0, 1)));  // "A"
    EXPECT_EQ(it, it_cpy);

    // extend_right(empty range)
    it_cpy = it;
    EXPECT_TRUE(it.extend_right(this->empty_text));
    EXPECT_EQ(it, it_cpy);
}

TYPED_TEST_P(fm_index_cursor_collection_test, extend_right_range_empty_text)
{
    typename TypeParam::index_type fm{this->text_col3}; // {"ACGACG", "", "", "TGCGATCGA"}

    // successful extend_right(range)
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text1, 1, 3)));   // "CG"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,1}, {0,4}, {3,2}, {3,6}}));
    EXPECT_EQ(it.query_length(), 2u);
    EXPECT_EQ(it.count(), 4u);

    EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text1, 0, 1)));   // "A"
    EXPECT_EQ(it.locate(), (std::vector<std::pair<uint64_t, uint64_t>>{{0,1}, {3,2}, {3,6}}));
    EXPECT_EQ(it.query_length(), 3u);
    EXPECT_EQ(it.count(), 3u);

    // unsuccessful extend_right(range), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.extend_right(seqan3::views::slice(this->text1, 0, 1)));  // "A"
    EXPECT_EQ(it, it_cpy);

    // extend_right(empty range)
    it_cpy = it;
    EXPECT_TRUE(it.extend_right(this->empty_text));
    EXPECT_EQ(it, it_cpy);
}

// TODO: doesn't work with the current structure of typed tests
// TYPED_TEST_P(fm_index_cursor_collection_test, extend_right_convertible_range)
// {
//     const std::vector<seqan3::dna4_vector> text{"ANGACGNN"_dna5};
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
    typename TypeParam::index_type fm{this->text_col2}; // {"ACGACG", "TGCGATCGA"}

    // successful extend_right(char)
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text1, 0, 1)));   // "A"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,0}, {0,3}, {1,4}, {1,8}}));
    EXPECT_EQ(it.query_length(), 1u);

    EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text1, 1, 2)));   // "C"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,0}, {0,3}}));
    EXPECT_EQ(it.query_length(), 2u);

    // unsuccessful extend_right(char), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.extend_right(seqan3::views::slice(this->text1, 1, 2)));  // "C"
    EXPECT_EQ(it, it_cpy);
}

// TODO: doesn't work with the current structure of typed tests
// TYPED_TEST_P(fm_index_cursor_collection_test, extend_right_convertible_char)
// {
//     const std::vector<seqan3::dna4_vector> text{"ANGACGNN"_dna5};
//     typename TypeParam::index_type fm{text};
//
//     // successful extend_right(char) using a different alphabet
//     TypeParam it(fm);
//     EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text1, 0, 1)));    // "A"
//     EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{0, 3}));
//     EXPECT_EQ(it.query_length(), 1);
// }

TYPED_TEST_P(fm_index_cursor_collection_test, extend_right_range_and_cycle)
{
    typename TypeParam::index_type fm{this->text_col4}; // {"ACGAACGC", "TACGATCGA"}

    // successful extend_right() and cycle_back()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text1, 0, 4)));   // "ACGA"
    EXPECT_EQ(it.locate(), (std::vector<std::pair<uint64_t, uint64_t>>{{0,0}, {1,1}}));
    EXPECT_EQ(it.query_length(), 4u);

    EXPECT_TRUE(it.cycle_back());
    EXPECT_EQ(it.locate(), (std::vector<std::pair<uint64_t, uint64_t>>{{0,4}}));
    EXPECT_EQ(it.query_length(), 4u);
}

TYPED_TEST_P(fm_index_cursor_collection_test, extend_right_char_and_cycle)
{
    typename TypeParam::index_type fm{this->text_col5}; // {"ACGAACGC"_dna4, "TGCGATCGA"_dna4};

    // successful extend_right() and cycle_back()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text1, 0, 1)));   // "A"
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,0}, {0,3}, {0,4},
                                                                                         {1,4}, {1,8}}));
    EXPECT_EQ(it.query_length(), 1u);

    EXPECT_TRUE(it.cycle_back());
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,1}, {0,5}, {0,7},
                                                                                         {1,2}, {1,6}}));
    EXPECT_EQ(it.query_length(), 1u);
}

TYPED_TEST_P(fm_index_cursor_collection_test, extend_right_and_cycle)
{
    typename TypeParam::index_type fm{this->text_col2}; // {"ACGACG", "TGCGATCGA"}

    // successful extend_right() and cycle_back()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right());
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,0}, {0,3}, {1,4}, {1,8}}));
    EXPECT_EQ(it.query_length(), 1u);

    EXPECT_TRUE(it.cycle_back());
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,1}, {0,4}, {1,2}, {1,6}}));
    EXPECT_EQ(it.query_length(), 1u);

    EXPECT_TRUE(it.extend_right());
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<std::pair<uint64_t, uint64_t>>{{0,1}, {0,4}, {1,2}, {1,6}}));
    EXPECT_EQ(it.query_length(), 2u);

    // unsuccessful cycle_back(), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.cycle_back());
    EXPECT_EQ(it, it_cpy);

    // unsuccessful extend_right(), it remains untouched
    it = TypeParam(fm);
    EXPECT_TRUE(it.extend_right(this->text1));   // "ACGACG"
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
    typename TypeParam::index_type fm{this->text_col2}; // {"ACGACG", "TGCGATCGA"}

    // query()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text1, 0, 3)));      // "ACG"
    EXPECT_TRUE(std::ranges::equal(it.path_label(this->text_col2),
                                   seqan3::views::slice(this->text1, 0, 3)));   // "ACG"
}

TYPED_TEST_P(fm_index_cursor_collection_test, last_rank)
{
    typename TypeParam::index_type fm{this->text_col2}; // {"ACGACG", "TGCGATCGA"}

    // last_rank()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right(seqan3::views::slice(this->text1, 0, 3)));  // "ACG"
    bool a = it.last_rank() == seqan3::to_rank(this->text1[2]);             // 'G'
    EXPECT_TRUE(a);
}

TYPED_TEST_P(fm_index_cursor_collection_test, incomplete_alphabet)
{
    using alphabet_type = typename TestFixture::alphabet_type;

    // search a char that does not occur in the text (higher rank than largest char occurring in text)
    {
        typename TypeParam::index_type fm{this->text_col1}; // {"ACGACG", "ACGACG"}
        TypeParam it = TypeParam(fm);
        EXPECT_FALSE(it.extend_right(seqan3::assign_char_to('T', alphabet_type{}))); // 'T'
        EXPECT_EQ(it, TypeParam(fm));
    }

    // search a char that does not occur in the text (smaller rank than smallest char occurring in text)
    {
        typename TypeParam::index_type fm{this->text_col6}; // {"CGTCGT", "CGTCGT"}
        TypeParam it = TypeParam(fm);
        EXPECT_FALSE(it.extend_right(seqan3::assign_rank_to(0, alphabet_type{})));   // 'A' or letter with smallest rank
        EXPECT_EQ(it, TypeParam(fm));
    }

    // search a char that does not occur in the text
    // (some rank that is neither the smallest nor the highest occurring in text)
    {
        typename TypeParam::index_type fm{this->text_col7}; // {"ATATAT", "ATATAT"}

        TypeParam it = TypeParam(fm);
        // get rank which is neither the smallest nor the highest:
        uint8_t middle_rank = std::round((seqan3::to_rank(this->text4[1]) +                     // 'T' and
                                          seqan3::to_rank(this->text4[0]))/2);                  // 'A'
        EXPECT_FALSE(it.extend_right(seqan3::assign_rank_to(middle_rank, alphabet_type{})));    // 'C'
        EXPECT_FALSE(it.extend_right(this->text1[2]));                                          // 'G'
        EXPECT_FALSE(it.extend_right(seqan3::views::slice(this->text7, 0, 4)));                 // "ACGT"
        EXPECT_FALSE(it.extend_right(seqan3::views::slice(this->text1, 2, 3)));                 // "G"
        EXPECT_EQ(it, TypeParam(fm));

        EXPECT_TRUE(it.extend_right(this->text4[0]));                                           // 'A'
        EXPECT_TRUE(it.cycle_back());
        EXPECT_TRUE(std::ranges::equal(it.path_label(this->text_col7),
                                       seqan3::views::slice(this->text4, 1, 2)));               // "T"
    }
}

TYPED_TEST_P(fm_index_cursor_collection_test, lazy_locate)
{
    typename TypeParam::index_type fm{this->text_col8}; // {"ACGTACGT", "TGCGATACGA"}

    TypeParam it = TypeParam(fm);
    it.extend_right(seqan3::views::slice(this->text1, 0, 3));    // "ACG"

    EXPECT_TRUE(std::ranges::equal(it.locate(), it.lazy_locate()));
}

TYPED_TEST_P(fm_index_cursor_collection_test, extend_const_char_pointer)
{
    using alphabet_type = typename TestFixture::alphabet_type;

    // test case for https://github.com/seqan/seqan3/issues/1473
    if constexpr (std::is_same<alphabet_type, char>::value)
    {
        typename TypeParam::index_type fm{this->text_col1}; // {"ACGACG", "ACGACG"}
        char const * cg = "CG";

        // extend_right()
        TypeParam it1 = TypeParam(fm);
        TypeParam it2 = TypeParam(fm);

        it1.extend_right(cg);
        it2.extend_right(seqan3::views::slice(this->text1, 1, 3));      // "CG"

        EXPECT_TRUE(std::ranges::equal(it1.locate(), it2.locate()));    // [(0,1),(0,4),(1,4),(1,1)]
    }
}

TYPED_TEST_P(fm_index_cursor_collection_test, concept_check)
{
    EXPECT_TRUE(seqan3::fm_index_cursor_specialisation<TypeParam>);
}

REGISTER_TYPED_TEST_SUITE_P(fm_index_cursor_collection_test, ctr, begin, extend_right_range,
                            extend_right_range_empty_text, extend_right_char, extend_right_range_and_cycle,
                            extend_right_char_and_cycle, extend_right_and_cycle, query, last_rank, incomplete_alphabet,
                            lazy_locate, extend_const_char_pointer, concept_check);
