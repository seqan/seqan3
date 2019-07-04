// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <type_traits>

#include "helper.hpp"

#include <seqan3/search/algorithm/all.hpp>

#include <gtest/gtest.h>

using namespace seqan3;
using namespace seqan3::search_cfg;
using namespace std::string_literals;

template <typename T>
class search_test : public ::testing::Test
{
public:
    std::vector<std::vector<dna4>> text{"ACGTACGTACGT"_dna4, "ACGTACGTACGT"_dna4};
    T index{text};
};

template <typename T>
class search_string_test : public ::testing::Test
{
public:
    std::vector<std::string> text{"Garfield the fat cat.", "Yet another text at position 1."};
    T index{text};
};

using fm_index_types        = ::testing::Types<fm_index<text_layout::collection>, bi_fm_index<text_layout::collection>>;
using fm_index_string_types = ::testing::Types<fm_index<text_layout::collection>, bi_fm_index<text_layout::collection>>;

TYPED_TEST_CASE(search_test, fm_index_types);
TYPED_TEST_CASE(search_string_test, fm_index_string_types);

TYPED_TEST(search_test, error_free)
{
    using result_t = std::pair<typename TypeParam::size_type, typename TypeParam::size_type>;
    using hits_result_t = std::vector<result_t>;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                             {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search with empty cfg
        configuration const cfg;
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                                  {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using empty max_total_error
        configuration const cfg = max_error{};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                                  {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using short version of max_total_error
        configuration const cfg = max_error{total{0}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                                  {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using max_total_error
        configuration const cfg = max_error{total{0}, substitution{0}, insertion{0}, deletion{0}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                                  {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using empty max_total_error_rate
        configuration const cfg = max_error_rate{};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                                  {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using short version of max_total_error_rate
        configuration const cfg = max_error_rate{total{.0}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                                  {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using max_total_error_rate
        configuration const cfg = max_error_rate{total{.0}, substitution{.0}, insertion{.0}, deletion{.0}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                                  {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }
}

TYPED_TEST(search_test, multiple_queries)
{
    using result_t = std::vector<std::pair<typename TypeParam::size_type, typename TypeParam::size_type>>;
    using hits_result_t = std::vector<result_t>;
    std::vector<std::vector<dna4>> const queries{{"GG"_dna4, "ACGTACGTACGT"_dna4, "ACGTA"_dna4}};

    configuration const cfg = max_error_rate{total{.0}, substitution{.0}, insertion{.0}, deletion{.0}};
    EXPECT_EQ(uniquify(search(queries, this->index, cfg)), (hits_result_t{{},
                                                                          {{0, 0}, {1, 0}},
                                                                          {{0, 0}, {0, 4}, {1, 0}, {1, 4}}}));
}

TYPED_TEST(search_string_test, error_free_string)
{
    using result_t = std::pair<typename TypeParam::size_type, typename TypeParam::size_type>;
    using hits_result_t = std::vector<result_t>;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search("at"s, this->index)), (hits_result_t{{0, 14}, {0, 18}, {1, 17}}));
        EXPECT_EQ(uniquify(search("Jon"s, this->index)), (hits_result_t{}));
    }
}

TYPED_TEST(search_string_test, error_free_raw)
{
    using result_t = std::pair<typename TypeParam::size_type, typename TypeParam::size_type>;
    using hits_result_t = std::vector<result_t>;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search("at", this->index)), (hits_result_t{{0, 14}, {0, 18}, {1, 17}}));
        EXPECT_EQ(uniquify(search("Jon", this->index)), (hits_result_t{}));
    }
}

TYPED_TEST(search_string_test, multiple_queries_string)
{
    using result_t = std::vector<std::pair<typename TypeParam::size_type, typename TypeParam::size_type>>;
    using hits_result_t = std::vector<result_t>;

    std::vector<std::string> const queries{"at", "Jon"};

    EXPECT_EQ(uniquify(search(queries, this->index)), (hits_result_t{{{0, 14}, {0, 18}, {1, 17}},
                                                                     {}})); // 3 and 0 hits
}

TYPED_TEST(search_string_test, multiple_queries_raw)
{
    using result_t = std::vector<std::pair<typename TypeParam::size_type, typename TypeParam::size_type>>;
    using hits_result_t = std::vector<result_t>;

    EXPECT_EQ(uniquify(search({"at", "Jon"}, this->index)), (hits_result_t{{{0, 14}, {0, 18}, {1, 17}},
                                                                           {}})); // 3 and 0 hits
}
