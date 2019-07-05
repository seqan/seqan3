// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <sstream>
#include <stdexcept>

#include <seqan3/alignment/multiple/library.hpp>
#include <seqan3/std/iterator>

using namespace seqan3::detail;

template <typename T>
struct library_test : public ::testing::Test
{
    msa_library<T> init_library_test()
    {
        msa_library<T> lib{};
        lib.insert({1, 2}, {5, 7}, 3);
        lib.insert({1, 2}, {5, 8}, 3);
        lib.insert({1, 3}, {5, 8}, 3);
        return lib;
    }
};

using test_types = ::testing::Types<int, double, size_t>;
TYPED_TEST_CASE(library_test, test_types);

TYPED_TEST(library_test, construction)
{
    EXPECT_TRUE(std::is_default_constructible<msa_library<TypeParam>>::value);
    EXPECT_TRUE(std::is_copy_constructible<msa_library<TypeParam>>::value);
    EXPECT_TRUE(std::is_copy_assignable<msa_library<TypeParam>>::value);
    EXPECT_TRUE(std::is_move_constructible<msa_library<TypeParam>>::value);
    EXPECT_TRUE(std::is_move_assignable<msa_library<TypeParam>>::value);
    EXPECT_TRUE(std::is_destructible<msa_library<TypeParam>>::value);
}

TYPED_TEST(library_test, insert_score_entry)
{
    msa_library<TypeParam> lib{};
    EXPECT_TRUE(lib.insert({1, 2}, {5, 7}, 3));
    EXPECT_TRUE(lib.insert({1, 2}, {5, 8}, 3));
    EXPECT_TRUE(lib.insert({1, 3}, {5, 8}, 3));
    EXPECT_FALSE(lib.insert({1, 2}, {5, 8}, 3));
    EXPECT_FALSE(lib.insert({1, 2}, {5, 8}, 1));
    EXPECT_FALSE(lib.insert({1, 2}, {5, 7}, 3));
    EXPECT_FALSE(lib.insert({1, 3}, {5, 8}, 3));
}

TYPED_TEST(library_test, add_score_entry)
{
    auto lib = this->init_library_test();

    // increase or decrease existing values
    auto elem = lib[{1, 2, 5, 7}]; // retrieve iterator to the element
    lib.add({1, 2}, {5, 7}, 10);
    EXPECT_EQ(elem.score(), static_cast<TypeParam>(13)); // 3 + 10
    lib.add({1, 2}, {5, 7}, -5);
    EXPECT_EQ(elem.score(), static_cast<TypeParam>(8));  // 13 - 5

    // add to non-existing entries
    EXPECT_EQ((lib[{1, 2, 4, 5}]), lib.end()); // does not exist
    lib.add({1, 2}, {4, 5}, 10);
    EXPECT_NE((lib[{1, 2, 4, 5}]), lib.end()); // exists now with score 10
    EXPECT_EQ((lib[{1, 2, 4, 5}].score()), static_cast<TypeParam>(10));
}

TYPED_TEST(library_test, member_access)
{
    auto lib = this->init_library_test();

    auto elem = lib[{1, 2, 5, 7}];
    EXPECT_NE(elem, lib.end()); // element exists
    EXPECT_EQ(elem.score(), static_cast<TypeParam>(3));
    EXPECT_EQ(elem.seq_pair(), std::make_pair(1ul, 2ul));
    EXPECT_EQ(elem.pos_pair(), std::make_pair(5ul, 7ul));
}

TYPED_TEST(library_test, alignment_access)
{
    auto lib = this->init_library_test();

    // element exists... check the values
    EXPECT_TRUE((lib[{1, 2}]));
    auto map_pos_score = lib[{1, 2}].value();
    EXPECT_EQ(map_pos_score.size(), 2ul);

    // examine the map of position pairs and score
    auto map_it = map_pos_score.begin();
    EXPECT_EQ(map_it->first, std::make_pair(5ul, 7ul));
    EXPECT_EQ(map_it->second, static_cast<TypeParam>(3));
    ++map_it;
    EXPECT_EQ(map_it->first, std::make_pair(5ul, 8ul));
    EXPECT_EQ(map_it->second, static_cast<TypeParam>(3));
    EXPECT_EQ(++map_it, map_pos_score.end());

    // element does not exist
    EXPECT_FALSE((lib[{1, 4}]));
    EXPECT_THROW((lib[{1, 4}].value()), std::bad_optional_access);
}

TYPED_TEST(library_test, empty_iterator)
{
    msa_library<TypeParam> lib{};
    auto empty_it = lib.begin();
    EXPECT_EQ(empty_it, lib.end()); // is at end
}

// TODO: use the iterator test template when available
TYPED_TEST(library_test, iterator)
{
    auto lib = this->init_library_test();

    // retrieve begin iterator
    auto it = lib.begin();
    EXPECT_TRUE(std::BidirectionalIterator<decltype(it)>);

    // increment
    ++it;

    // verify the values
    auto && [seq, pos, score] = *it;
    EXPECT_EQ(seq, std::make_pair(1ul, 2ul));
    EXPECT_EQ(pos, std::make_pair(5ul, 8ul));
    EXPECT_EQ(score, static_cast<TypeParam>(3));

    // assign to score and verify
    score = 20;
    EXPECT_EQ((lib[{1, 2, 5, 8}].score()), static_cast<TypeParam>(20));

    // increment with switch to next sequence pair
    it++;
    EXPECT_EQ(it.seq_pair(), std::make_pair(1ul, 3ul));
    EXPECT_EQ(it.pos_pair(), std::make_pair(5ul, 8ul));
    EXPECT_EQ(it.score(), static_cast<TypeParam>(3));

    // reach the end
    ++it;
    EXPECT_EQ(it, lib.end());

    // decrement until begin
    it--;
    --it;
    --it;
    EXPECT_EQ(it, lib.begin());
}

TYPED_TEST(library_test, stream)
{
    auto lib = this->init_library_test();

    std::ostringstream stream;
    stream << lib;
    EXPECT_EQ(stream.str(), "# 1 2\n"
                            "5 7 3\n"
                            "5 8 3\n"
                            "# 1 3\n"
                            "5 8 3\n");
}

TYPED_TEST(library_test, lib_format)
{
    auto lib = this->init_library_test();
    std::vector<std::string> ids{"id_one", "id_two"};
    std::vector<std::string> seqs{"GCGCUUAGCAA", "UUGCUCGAAGCC"};

    std::ostringstream stream{};
    stream << std::make_tuple(lib, ids, seqs);
    EXPECT_EQ(stream.str(), "! T-COFFEE_LIB_FORMAT_01\n"
                            "2\n"
                            "id_one 11 GCGCUUAGCAA\n"
                            "id_two 12 UUGCUCGAAGCC\n"
                            "# 1 2\n"
                            "5 7 3\n"
                            "5 8 3\n"
                            "# 1 3\n"
                            "5 8 3\n"
                            "! SEQ_1_TO_N\n");
}
