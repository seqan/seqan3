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
    std::vector<dna4> text{"ACGTACGTACGT"_dna4};
    T index{text};
};

template <typename T>
class search_string_test : public ::testing::Test
{
public:
    std::string text{"Garfield the fat cat."};
    T index{text};
};

using fm_index_types        = ::testing::Types<fm_index<text_layout::single>, bi_fm_index<text_layout::single>>;
using fm_index_string_types = ::testing::Types<fm_index<text_layout::single>, bi_fm_index<text_layout::single>>;

TYPED_TEST_CASE(search_test, fm_index_types);
TYPED_TEST_CASE(search_string_test, fm_index_string_types);

TYPED_TEST(search_test, error_free)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search with empty cfg
        configuration const cfg;
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using empty max_total_error
        configuration const cfg = max_error{};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using short version of max_total_error
        configuration const cfg = max_error{total{0}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using max_total_error
        configuration const cfg = max_error{total{0}, substitution{0}, insertion{0}, deletion{0}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using empty max_total_error_rate
        configuration const cfg = max_error_rate{};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using short version of max_total_error_rate
        configuration const cfg = max_error_rate{total{.0}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using max_total_error_rate
        configuration const cfg = max_error_rate{total{.0}, substitution{.0}, insertion{.0}, deletion{.0}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }
}

TYPED_TEST(search_test, multiple_queries)
{
    using hits_result_t = std::vector<std::vector<typename TypeParam::size_type>>;
    std::vector<std::vector<dna4>> const queries{{"GG"_dna4, "ACGTACGTACGT"_dna4, "ACGTA"_dna4}};

    configuration const cfg = max_error_rate{total{.0}, substitution{.0}, insertion{.0}, deletion{.0}};
    EXPECT_EQ(uniquify(search(queries, this->index, cfg)), (hits_result_t{{}, {0}, {0, 4}})); // 0, 1 and 2 hits
}

TYPED_TEST(search_test, invalid_error_configuration)
{
    configuration const cfg = max_error{total{0}, substitution{1}};
    EXPECT_THROW(search("A"_dna4, this->index, cfg), std::invalid_argument);
}

TYPED_TEST(search_test, error_substitution)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error_rate{total{.25}, substitution{.25}};

        EXPECT_EQ(uniquify(search("ACGT"_dna4    , this->index, cfg)), (hits_result_t{0, 4, 8})); // exact match
        EXPECT_EQ(uniquify(search("CGG"_dna4     , this->index, cfg)), (hits_result_t{}));        // not enough mismatches
        EXPECT_EQ(uniquify(search("CGTC"_dna4    , this->index, cfg)), (hits_result_t{1, 5}));    // 1 mismatch
        EXPECT_EQ(uniquify(search("ACGGACG"_dna4 , this->index, cfg)), (hits_result_t{0, 4}));    // 1 mismatch
        EXPECT_EQ(uniquify(search("ACGGACGG"_dna4, this->index, cfg)), (hits_result_t{0, 4}));    // 2 mismatches
    }

    {
        configuration const cfg = max_error_rate{total{.25}, substitution{.25}, insertion{.0}, deletion{.0}};

        EXPECT_EQ(uniquify(search("ACGT"_dna4    , this->index, cfg)), (hits_result_t{0, 4, 8})); // exact match
        EXPECT_EQ(uniquify(search("CGG"_dna4     , this->index, cfg)), (hits_result_t{}));        // not enough mismatches
        EXPECT_EQ(uniquify(search("CGTC"_dna4    , this->index, cfg)), (hits_result_t{1, 5}));    // 1 mismatch
        EXPECT_EQ(uniquify(search("ACGGACG"_dna4 , this->index, cfg)), (hits_result_t{0, 4}));    // 1 mismatch
        EXPECT_EQ(uniquify(search("ACGGACGG"_dna4, this->index, cfg)), (hits_result_t{0, 4}));    // 2 mismatches
    }

    {
        configuration const cfg = max_error{total{1}, substitution{1}};

        EXPECT_EQ(uniquify(search("ACGT"_dna4    , this->index, cfg)), (hits_result_t{0, 4, 8})); // exact match
        EXPECT_EQ(uniquify(search("CGTTT"_dna4   , this->index, cfg)), (hits_result_t{}));        // not enough mismatches
        EXPECT_EQ(uniquify(search("CGG"_dna4     , this->index, cfg)), (hits_result_t{1, 5, 9})); // 1 mismatch
        EXPECT_EQ(uniquify(search("ACGGACG"_dna4 , this->index, cfg)), (hits_result_t{0, 4}));    // 1 mismatch
        EXPECT_EQ(uniquify(search("CGTCCGTA"_dna4, this->index, cfg)), (hits_result_t{1}));       // 1 mismatch
    }

    {
        configuration const cfg = max_error{total{1}, substitution{1}, insertion{0}, deletion{0}};

        EXPECT_EQ(uniquify(search("ACGT"_dna4    , this->index, cfg)), (hits_result_t{0, 4, 8})); // exact match
        EXPECT_EQ(uniquify(search("CGTTT"_dna4   , this->index, cfg)), (hits_result_t{}));        // not enough mismatches
        EXPECT_EQ(uniquify(search("CGG"_dna4     , this->index, cfg)), (hits_result_t{1, 5, 9})); // 1 mismatch
        EXPECT_EQ(uniquify(search("ACGGACG"_dna4 , this->index, cfg)), (hits_result_t{0, 4}));    // 1 mismatch
        EXPECT_EQ(uniquify(search("CGTCCGTA"_dna4, this->index, cfg)), (hits_result_t{1}));       // 1 mismatch
    }
}

TYPED_TEST(search_test, error_configuration_types)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        uint8_t s = 1, t = 1;
        configuration const cfg = max_error{total{t}, substitution{s}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
    }

    {
        configuration const cfg = max_error{substitution{1}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
    }
}

TYPED_TEST(search_test, error_insertion)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error_rate{total{.25}, insertion{.25}};

        // exact match and insertion at the beginning of the query
        EXPECT_EQ(uniquify(search("ACGT"_dna4    , this->index, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
        // 1 insertion
        EXPECT_EQ(uniquify(search("CCGT"_dna4    , this->index, cfg)), (hits_result_t{1, 5, 9}));
        // 2 insertions
        EXPECT_EQ(uniquify(search("ACCGGTAC"_dna4, this->index, cfg)), (hits_result_t{0, 4}));
        // 2 insertions necessary, only 1 allowed
        EXPECT_EQ(uniquify(search("ACCGG"_dna4   , this->index, cfg)), (hits_result_t{}));
        // deletion necessary, not allowed
        EXPECT_EQ(uniquify(search("ACTACGT"_dna4 , this->index, cfg)), (hits_result_t{}));
    }

    {
        configuration const cfg = max_error{total{1}, insertion{1}};

        // exact match and insertion at the beginning of the query
        EXPECT_EQ(uniquify(search("ACGT"_dna4    , this->index, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
        // 1 insertion
        EXPECT_EQ(uniquify(search("CCGT"_dna4    , this->index, cfg)), (hits_result_t{1, 5, 9}));
        // 2 insertions necessary, only 1 allowed
        EXPECT_EQ(uniquify(search("ACCGGTAC"_dna4, this->index, cfg)), (hits_result_t{}));
        // deletion necessary, not allowed
        EXPECT_EQ(uniquify(search("ACTACGT"_dna4 , this->index, cfg)), (hits_result_t{}));
    }
}

TYPED_TEST(search_test, error_deletion)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error_rate{total{.25}, deletion{.25}};

        // exact match, no deletion
        EXPECT_EQ(uniquify(search("ACGT"_dna4    , this->index, cfg)), (hits_result_t{0, 4, 8}));
        // not enough max errors
        EXPECT_EQ(uniquify(search("AGT"_dna4     , this->index, cfg)), (hits_result_t{}));
        // one deletion (C)
        EXPECT_EQ(uniquify(search("AGTA"_dna4    , this->index, cfg)), (hits_result_t{0, 4}));
        // two deletion (C)
        EXPECT_EQ(uniquify(search("AGTAGTAC"_dna4, this->index, cfg)), (hits_result_t{0}));
        // no deletion at beginning. 0 and 4 cannot be reported
        EXPECT_EQ(uniquify(search("CGTACGT"_dna4 , this->index, cfg)), (hits_result_t{1, 5}));
    }

    {
        configuration const cfg = max_error{total{1}, deletion{1}};

        // exact match, no deletion
        EXPECT_EQ(uniquify(search("ACGT"_dna4    , this->index, cfg)), (hits_result_t{0, 4, 8}));
        // one deletion (C)
        EXPECT_EQ(uniquify(search("AGTA"_dna4    , this->index, cfg)), (hits_result_t{0, 4}));
        // 2 deletions necessary, only 1 allowed
        EXPECT_EQ(uniquify(search("AGTAGTAC"_dna4, this->index, cfg)), (hits_result_t{}));
        // no deletion at beginning. 0 and 4 cannot be reported
        EXPECT_EQ(uniquify(search("CGTACGT"_dna4 , this->index, cfg)), (hits_result_t{1, 5}));
    }
}

TYPED_TEST(search_test, error_levenshtein)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error{total{1}};
        EXPECT_EQ(uniquify(search("CCGT"_dna4, this->index, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }

    {
        configuration const cfg = max_error{total{2}};
        EXPECT_EQ(uniquify(search("CCGT"_dna4, this->index, cfg)), (hits_result_t{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
    }
}

TYPED_TEST(search_test, search_strategy_all)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error{total{1}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }

    {
        configuration const cfg = max_error{total{1}} | mode{all};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }
}

TYPED_TEST(search_test, search_strategy_best)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error{total{1}} | mode{best};

        hits_result_t possible_hits{0, 4, 8}; // any of 0, 4, 8 ... 1, 5, 9 are not best hits
        hits_result_t result = search("ACGT"_dna4, this->index, cfg);
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());

        EXPECT_EQ(search("AAAA"_dna4, this->index, cfg), (hits_result_t{})); // no hit
    }
}

TYPED_TEST(search_test, search_strategy_all_best)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error{total{1}} | mode{all_best};

        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8})); // 1, 5, 9 are not best hits

        EXPECT_EQ(search("AAAA"_dna4, this->index, cfg), (hits_result_t{})); // no hit
    }
}

TYPED_TEST(search_test, search_strategy_strata)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        configuration const cfg = max_error{total{1}} | mode{strata{0}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
    }

    {
        configuration const cfg = max_error{total{1}} | mode{strata{1}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }

    {
        configuration const cfg = max_error{total{1}} | mode{strata{1}};
        EXPECT_EQ(search("AAAA"_dna4, this->index, cfg), (hits_result_t{})); // no hit
    }

    // {
    //     // best hit ACGT with 1 error, i.e. 1+1
    //     configuration const cfg = max_total_error{1} | strategy_strata{1};
    //     EXPECT_EQ(uniquify(search(this->index, "CCGT"_dna4, cfg)), (hits_result_t{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
    // }

    // {
    //     // best hit ACGT with 1 error, i.e. 1+1
    //     configuration const cfg = max_total_error{1} | strategy_strata{1};
    //     EXPECT_EQ(uniquify(search(this->index, "CCGT"_dna4, cfg)), (hits_result_t{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
    // }
}

TYPED_TEST(search_string_test, error_free_string)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search("at"s, this->index)), (hits_result_t{14, 18}));
        EXPECT_EQ(uniquify(search("Jon"s, this->index)), (hits_result_t{}));
    }
}

TYPED_TEST(search_string_test, error_free_raw)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search("at", this->index)), (hits_result_t{14, 18}));
        EXPECT_EQ(uniquify(search("Jon", this->index)), (hits_result_t{}));
    }
}

TYPED_TEST(search_string_test, multiple_queries_string)
{
    using hits_result_t = std::vector<std::vector<typename TypeParam::size_type>>;
    std::vector<std::string> const queries{"at", "Jon"};

    EXPECT_EQ(uniquify(search(queries, this->index)), (hits_result_t{{14, 18}, {}})); // 2 and 0 hits
}

TYPED_TEST(search_string_test, multiple_queries_raw)
{
    using hits_result_t = std::vector<std::vector<typename TypeParam::size_type>>;

    EXPECT_EQ(uniquify(search({"at", "Jon"}, this->index)), (hits_result_t{{14, 18}, {}})); // 2 and 0 hits
}

// TYPED_TEST(search_test, return_iterator_index)
// {
// }
//
// TYPED_TEST(search_test, on_hit)
// {
// }
