// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <type_traits>

#include "helper.hpp"

#include <seqan3/search/algorithm/all.hpp>

#include <gtest/gtest.h>

using namespace seqan3;
using namespace seqan3::search_cfg;

template <typename T>
class search_test : public ::testing::Test
{
public:
    std::vector<std::vector<dna4>> text{"ACGTACGTACGT"_dna4, "ACGTACGTACGT"_dna4};
    T index{text};
};

using fm_index_types = ::testing::Types<fm_index<std::vector<std::vector<dna4>>>,
                                        bi_fm_index<std::vector<std::vector<dna4>>>>;

TYPED_TEST_CASE(search_test, fm_index_types);

TYPED_TEST(search_test, error_free)
{
    using result_t = std::pair<typename TypeParam::size_type, typename TypeParam::size_type>;
    using hits_result_t = std::vector<result_t>;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                             {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search with empty cfg
        configuration const cfg;
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                                  {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using empty max_total_error
        configuration const cfg = max_error{};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                                  {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using short version of max_total_error
        configuration const cfg = max_error{total{0}};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                                  {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using max_total_error
        configuration const cfg = max_error{total{0}, substitution{0}, insertion{0}, deletion{0}};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                                  {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using empty max_total_error_rate
        configuration const cfg = max_error_rate{};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                                  {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using short version of max_total_error_rate
        configuration const cfg = max_error_rate{total{.0}};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                                  {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using max_total_error_rate
        configuration const cfg = max_error_rate{total{.0}, substitution{.0}, insertion{.0}, deletion{.0}};
        EXPECT_EQ(uniquify(search(this->index, "ACGT"_dna4, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8},
                                                                                  {1, 0}, {1, 4}, {1, 8}}));
        EXPECT_EQ(uniquify(search(this->index, "ACGG"_dna4, cfg)), (hits_result_t{}));
    }
}

TYPED_TEST(search_test, multiple_queries)
{
    using result_t = std::vector<std::pair<typename TypeParam::size_type, typename TypeParam::size_type>>;
    using hits_result_t = std::vector<result_t>;
    std::vector<std::vector<dna4>> const queries{{"GG"_dna4, "ACGTACGTACGT"_dna4, "ACGTA"_dna4}};

    configuration const cfg = max_error_rate{total{.0}, substitution{.0}, insertion{.0}, deletion{.0}};
    EXPECT_EQ(uniquify(search(this->index, queries, cfg)), (hits_result_t{{},
                                                                          {{0, 0}, {1, 0}},
                                                                          {{0, 0}, {0, 4}, {1, 0}, {1, 4}}}));
}
