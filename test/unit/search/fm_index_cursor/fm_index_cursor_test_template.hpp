/// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include "../helper.hpp"

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/std/algorithm>

#include <gtest/gtest.h>

using seqan3::operator""_dna4;

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
class fm_index_cursor_test : public ::testing::Test
{};

TYPED_TEST_SUITE_P(fm_index_cursor_test);

TYPED_TEST_P(fm_index_cursor_test, ctr)
{
    seqan3::dna4_vector text{"ACGACG"_dna4};
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

TYPED_TEST_P(fm_index_cursor_test, begin)
{
    seqan3::dna4_vector text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // begin
    TypeParam it(fm);
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{0, 1, 2, 3, 4, 5, 6}));// sentinel position included
    EXPECT_EQ(it.query_length(), 0u);
    EXPECT_EQ(it.count(), 7u);
}

TYPED_TEST_P(fm_index_cursor_test, extend_right_range)
{
    seqan3::dna4_vector text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right(range)
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right("CG"_dna4));
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{1, 4}));
    EXPECT_EQ(it.query_length(), 2u);
    EXPECT_EQ(it.count(), 2u);

    EXPECT_TRUE(it.extend_right("A"_dna4));
    EXPECT_EQ(it.locate(), (std::vector<uint64_t>{1}));
    EXPECT_EQ(it.query_length(), 3u);
    EXPECT_EQ(it.count(), 1u);

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
// TYPED_TEST_P(fm_index_cursor_test, extend_right_convertible_range)
// {
//     seqan3::dna4_vector text{"ANGACGNN"_dna5};
//     typename TypeParam::index_type fm{text};
//
//     // successful extend_right(range) using a different alphabet
//     TypeParam it(fm);
//     EXPECT_TRUE(it.extend_right("GA"_dna4));
//     EXPECT_EQ(it.locate(), (std::vector<uint64_t>{2}));
//     EXPECT_EQ(it.query_length(), 2);
// }

TYPED_TEST_P(fm_index_cursor_test, extend_right_char)
{
    seqan3::dna4_vector text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right(char)
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right('A'_dna4));
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{0, 3}));
    EXPECT_EQ(it.query_length(), 1u);

    EXPECT_TRUE(it.extend_right('C'_dna4));
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{0, 3}));
    EXPECT_EQ(it.query_length(), 2u);

    // unsuccessful extend_right(char), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.extend_right('C'_dna4));
    EXPECT_EQ(it, it_cpy);
}

// TODO: doesn't work with the current structure of typed tests
// TYPED_TEST_P(fm_index_cursor_test, extend_right_convertible_char)
// {
//     seqan3::dna4_vector text{"ANGACGNN"_dna5};
//     typename TypeParam::index_type fm{text};
//
//     // successful extend_right(char) using a different alphabet
//     TypeParam it(fm);
//     EXPECT_TRUE(it.extend_right('A'_dna4));
//     EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{0, 3}));
//     EXPECT_EQ(it.query_length(), 1);
// }

TYPED_TEST_P(fm_index_cursor_test, extend_right_range_and_cycle)
{
    seqan3::dna4_vector text{"ACGAACGC"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right() and cycle_back()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right("ACGA"_dna4));
    EXPECT_EQ(it.locate(), (std::vector<uint64_t>{0}));
    EXPECT_EQ(it.query_length(), 4u);

    EXPECT_TRUE(it.cycle_back());
    EXPECT_EQ(it.locate(), (std::vector<uint64_t>{4}));
    EXPECT_EQ(it.query_length(), 4u);
}

TYPED_TEST_P(fm_index_cursor_test, extend_right_char_and_cycle)
{
    seqan3::dna4_vector text{"ACGAACGC"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right() and cycle_back()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right('A'_dna4));
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{0, 3, 4}));
    EXPECT_EQ(it.query_length(), 1u);

    EXPECT_TRUE(it.cycle_back());
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{1, 5, 7}));
    EXPECT_EQ(it.query_length(), 1u);
}

TYPED_TEST_P(fm_index_cursor_test, extend_right_and_cycle)
{
    seqan3::dna4_vector text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // successful extend_right() and cycle_back()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right());
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{0, 3}));
    EXPECT_EQ(it.query_length(), 1u);

    EXPECT_TRUE(it.cycle_back());
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{1, 4}));
    EXPECT_EQ(it.query_length(), 1u);

    EXPECT_TRUE(it.extend_right());
    EXPECT_EQ(seqan3::uniquify(it.locate()), (std::vector<uint64_t>{1, 4}));
    EXPECT_EQ(it.query_length(), 2u);

    // unsuccessful cycle_back(), it remains untouched
    TypeParam it_cpy = it;
    EXPECT_FALSE(it.cycle_back());
    EXPECT_EQ(it, it_cpy);

    // unsuccessful extend_right(), it remains untouched
    it = TypeParam(fm);
    EXPECT_TRUE(it.extend_right("GACG"_dna4));
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

TYPED_TEST_P(fm_index_cursor_test, query)
{
    seqan3::dna4_vector text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // query()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right("ACG"_dna4));
    EXPECT_TRUE(std::ranges::equal(it.path_label(text), "ACG"_dna4));
}

TYPED_TEST_P(fm_index_cursor_test, last_rank)
{
    seqan3::dna4_vector text{"ACGACG"_dna4};
    typename TypeParam::index_type fm{text};

    // last_rank()
    TypeParam it(fm);
    EXPECT_TRUE(it.extend_right("ACG"_dna4));
    bool a = it.last_rank() == seqan3::to_rank('G'_dna4);
    EXPECT_TRUE(a);
}

TYPED_TEST_P(fm_index_cursor_test, incomplete_alphabet)
{
    // search a char that does not occur in the text (higher rank than largest char occurring in text)
    {
        seqan3::dna4_vector text{"ACGACG"_dna4};
        typename TypeParam::index_type fm{text};
        TypeParam it = TypeParam(fm);
        EXPECT_FALSE(it.extend_right('T'_dna4));
        EXPECT_EQ(it, TypeParam(fm));
    }

    // search a char that does not occur in the text (smaller rank than smallest char occurring in text)
    {
        seqan3::dna4_vector text{"CGTCGT"_dna4};
        typename TypeParam::index_type fm{text};
        TypeParam it = TypeParam(fm);
        EXPECT_FALSE(it.extend_right('A'_dna4));
        EXPECT_EQ(it, TypeParam(fm));
    }

    // search a char that does not occur in the text
    // (some rank that is neither the smallest nor the highest occurring in text)
    {
        seqan3::dna4_vector text{"ATATAT"_dna4};
        typename TypeParam::index_type fm{text};
        TypeParam it = TypeParam(fm);
        EXPECT_FALSE(it.extend_right('C'_dna4));
        EXPECT_FALSE(it.extend_right('G'_dna4));
        EXPECT_FALSE(it.extend_right("ACGT"_dna4));
        EXPECT_FALSE(it.extend_right("G"_dna4));
        EXPECT_EQ(it, TypeParam(fm));

        EXPECT_TRUE(it.extend_right('A'_dna4));
        EXPECT_TRUE(it.cycle_back());
        EXPECT_TRUE(std::ranges::equal(it.path_label(text), "T"_dna4));
    }
}

TYPED_TEST_P(fm_index_cursor_test, lazy_locate)
{
    seqan3::dna4_vector text{"ACGTACGT"_dna4};
    typename TypeParam::index_type fm{text};

    TypeParam it = TypeParam(fm);
    it.extend_right("ACG"_dna4);

    EXPECT_TRUE(std::ranges::equal(it.locate(), it.lazy_locate()));
}

TYPED_TEST_P(fm_index_cursor_test, concept_check)
{
    EXPECT_TRUE(seqan3::fm_index_cursor_specialisation<TypeParam>);
}

REGISTER_TYPED_TEST_SUITE_P(fm_index_cursor_test, ctr, begin, extend_right_range, extend_right_char,
                            extend_right_range_and_cycle, extend_right_char_and_cycle, extend_right_and_cycle, query,
                            last_rank, incomplete_alphabet, lazy_locate, concept_check);
