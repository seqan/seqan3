// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <type_traits>

#include <seqan3/search/all.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include "helper.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_phred42;

using namespace std::string_literals;

template <typename index_t>
class search_test : public ::testing::Test
{
public:
    using hits_result_t = std::vector<std::pair<size_t, typename index_t::size_type>>;

    seqan3::dna4_vector text{"ACGTACGTACGT"_dna4};
    index_t index{text};
};

template <typename index_t>
class search_string_test : public ::testing::Test
{
public:
    using hits_result_t = std::vector<std::pair<size_t, typename index_t::size_type>>;

    std::string text{"Garfield the fat cat."};
    index_t index{text};
};

using fm_index_types        = ::testing::Types<seqan3::fm_index<seqan3::dna4, seqan3::text_layout::single>,
                                               seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::single>>;
using fm_index_string_types = ::testing::Types<seqan3::fm_index<char, seqan3::text_layout::single>,
                                               seqan3::bi_fm_index<char, seqan3::text_layout::single>>;

TYPED_TEST_SUITE(search_test, fm_index_types, );
TYPED_TEST_SUITE(search_string_test, fm_index_string_types, );

TYPED_TEST(search_test, error_free)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index)), (hits_result_t{{0, 0}, {0, 4}, {0, 8}}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search with empty cfg
        seqan3::configuration const cfg;
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8}}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using empty max_total_error
        // default max_error{} sets all error to 0
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8}}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using empty max_total_error_rate
        // default max_error_rate{} sets all error to 0.0
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8}}));
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }
}

TYPED_TEST(search_test, convertible_query)
{
    using hits_result_t = typename TestFixture::hits_result_t;
    std::vector<seqan3::qualified<seqan3::dna4, seqan3::phred42>> query{{'A'_dna4, '!'_phred42},
                                                                        {'C'_dna4, '!'_phred42},
                                                                        {'G'_dna4, '!'_phred42},
                                                                        {'T'_dna4, '!'_phred42}};

    EXPECT_EQ(uniquify(search(query, this->index)), (hits_result_t{{0, 0}, {0, 4}, {0, 8}}));
}

TYPED_TEST(search_test, multiple_queries)
{
    using hits_result_t = typename TestFixture::hits_result_t;
    std::vector<std::vector<seqan3::dna4>> const queries{{"GG"_dna4, "ACGTACGTACGT"_dna4, "ACGTA"_dna4}};

    seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::total{.0},
                                                                         seqan3::search_cfg::substitution{.0},
                                                                         seqan3::search_cfg::insertion{.0},
                                                                         seqan3::search_cfg::deletion{.0}};
    EXPECT_EQ(uniquify(search(queries, this->index, cfg)), (hits_result_t{{1, 0}, {2, 0}, {2, 4}})); // 0, 1 and 2 hits
}

TYPED_TEST(search_test, invalid_error_configuration)
{
    seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{0},
                                                                            seqan3::search_cfg::substitution{1}};
    EXPECT_THROW(search("A"_dna4, this->index, cfg), std::invalid_argument);
}

TYPED_TEST(search_test, error_substitution)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::substitution{.25}};

        EXPECT_EQ(uniquify(search("ACGT"_dna4    , this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8}})); // exact match
        EXPECT_EQ(uniquify(search("CGG"_dna4     , this->index, cfg)), (hits_result_t{}));                       // not enough mismatches
        EXPECT_EQ(uniquify(search("CGTC"_dna4    , this->index, cfg)), (hits_result_t{{0, 1}, {0, 5}}));         // 1 mismatch
        EXPECT_EQ(uniquify(search("ACGGACG"_dna4 , this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}}));         // 1 mismatch
        EXPECT_EQ(uniquify(search("ACGGACGG"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}}));         // 2 mismatches
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::substitution{1}};

        EXPECT_EQ(uniquify(search("ACGT"_dna4    , this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8}})); // exact match
        EXPECT_EQ(uniquify(search("CGTTT"_dna4   , this->index, cfg)), (hits_result_t{}));                       // not enough mismatches
        EXPECT_EQ(uniquify(search("CGG"_dna4     , this->index, cfg)), (hits_result_t{{0, 1}, {0, 5}, {0, 9}})); // 1 mismatch
        EXPECT_EQ(uniquify(search("ACGGACG"_dna4 , this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}}));         // 1 mismatch
        EXPECT_EQ(uniquify(search("CGTCCGTA"_dna4, this->index, cfg)), (hits_result_t{{0, 1}}));                 // 1 mismatch
    }
}

TYPED_TEST(search_test, error_configuration_types)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::substitution{1}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8}}));
    }
}

TYPED_TEST(search_test, error_insertion)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::insertion{.25}};

        // exact match and insertion at the beginning of the query
        EXPECT_EQ(uniquify(search("ACGT"_dna4    , this->index, cfg)), (hits_result_t{{0, 0}, {0, 1}, {0, 4}, {0, 5},
                                                                                      {0, 8}, {0, 9}}));
        // 1 insertion
        EXPECT_EQ(uniquify(search("CCGT"_dna4    , this->index, cfg)), (hits_result_t{{0, 1}, {0, 5}, {0, 9}}));
        // 2 insertions
        EXPECT_EQ(uniquify(search("ACCGGTAC"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}}));
        // 2 insertions necessary, only 1 allowed
        EXPECT_EQ(seqan3::uniquify(search("ACCGG"_dna4,    this->index, cfg)), (hits_result_t{}));
        // deletion necessary, not allowed
        EXPECT_EQ(seqan3::uniquify(search("ACTACGT"_dna4,  this->index, cfg)), (hits_result_t{}));
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::insertion{1}};

        // exact match and insertion at the beginning of the query
        EXPECT_EQ(uniquify(search("ACGT"_dna4    , this->index, cfg)), (hits_result_t{{0, 0}, {0, 1}, {0, 4}, {0, 5},
                                                                                      {0, 8}, {0, 9}}));
        // 1 insertion
        EXPECT_EQ(uniquify(search("CCGT"_dna4    , this->index, cfg)), (hits_result_t{{0, 1}, {0, 5}, {0, 9}}));
        // 2 insertions necessary, only 1 allowed
        EXPECT_EQ(seqan3::uniquify(search("ACCGGTAC"_dna4, this->index, cfg)), (hits_result_t{}));
        // deletion necessary, not allowed
        EXPECT_EQ(seqan3::uniquify(search("ACTACGT"_dna4,  this->index, cfg)), (hits_result_t{}));
    }
}

TYPED_TEST(search_test, error_deletion)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::deletion{.25}};

        // exact match, no deletion
        EXPECT_EQ(uniquify(search("ACGT"_dna4    , this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8}}));
        // not enough max errors
        EXPECT_EQ(seqan3::uniquify(search("AGT"_dna4,      this->index, cfg)), (hits_result_t{}));
        // one deletion (C)
        EXPECT_EQ(uniquify(search("AGTA"_dna4    , this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}}));
        // two deletion (C)
        EXPECT_EQ(uniquify(search("AGTAGTAC"_dna4, this->index, cfg)), (hits_result_t{{0, 0}}));
        // no deletion at beginning. 0 and 4 cannot be reported
        EXPECT_EQ(uniquify(search("CGTACGT"_dna4 , this->index, cfg)), (hits_result_t{{0, 1}, {0, 5}}));
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::deletion{1}};

        // exact match, no deletion
        EXPECT_EQ(uniquify(search("ACGT"_dna4    , this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8}}));
        // one deletion (C)
        EXPECT_EQ(uniquify(search("AGTA"_dna4    , this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}}));
        // 2 deletions necessary, only 1 allowed
        EXPECT_EQ(seqan3::uniquify(search("AGTAGTAC"_dna4, this->index, cfg)), (hits_result_t{}));
        // no deletion at beginning. 0 and 4 cannot be reported
        EXPECT_EQ(uniquify(search("CGTACGT"_dna4 , this->index, cfg)), (hits_result_t{{0, 1}, {0, 5}}));
    }
}

TYPED_TEST(search_test, error_levenshtein)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}};
        EXPECT_EQ(uniquify(search("CCGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 1}, {0, 4}, {0, 5},
                                                                                  {0, 8}, {0, 9}}));
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{2}};
        EXPECT_EQ(uniquify(search("CCGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 1}, {0, 2}, {0, 3},
                                                                                  {0, 4}, {0, 5}, {0, 6}, {0, 7},
                                                                                  {0, 8}, {0, 9}, {0, 10}}));
    }
}

TYPED_TEST(search_test, error_indel_no_substitution)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    {
        // Match one mismatch with 1 insertion and deletion since mismatches are not allowed.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{2},
                                                                        seqan3::search_cfg::deletion{2},
                                                                        seqan3::search_cfg::insertion{2}};
        EXPECT_EQ(uniquify(search("GTACCTAC"_dna4, this->index, cfg)), (hits_result_t{{0, 2}}));
    }

    {
        // Enumerate a deletion and match one mismatch with 1 insertion and deletion since mismatches are not allowed.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{3},
                                                                        seqan3::search_cfg::deletion{3},
                                                                        seqan3::search_cfg::insertion{3}};
        EXPECT_EQ(uniquify(search("GTATCCTAC"_dna4, this->index, cfg)), (hits_result_t{{0, 2}}));
    }
}

TYPED_TEST(search_test, search_strategy_all)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 1}, {0, 4},
                                                                                  {0, 5}, {0, 8}, {0, 9}}));
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}} |
                                          seqan3::search_cfg::mode{seqan3::search_cfg::all};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 1}, {0, 4},
                                                                                   {0, 5}, {0, 8}, {0, 9}}));
    }
}

TYPED_TEST(search_test, search_strategy_best)
{
    using hits_result_t = typename TestFixture::hits_result_t;
    auto search_cfg_mode_best = seqan3::search_cfg::mode{seqan3::search_cfg::best};

    hits_result_t possible_hits{{0, 0}, {0, 4}, {0, 8}}; // any of 0, 4, 8 ... 1, 5, 9 are not best hits
    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}} |
                                          search_cfg_mode_best;

        hits_result_t result = search("ACGT"_dna4, this->index, cfg) | seqan3::views::to<std::vector>;
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());

        EXPECT_EQ(search("AAAA"_dna4, this->index, cfg) | seqan3::views::to<std::vector>, (hits_result_t{})); // no hit
    }

    { // Find best match with 1 insertion at the end.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::insertion{1}} |
                                          search_cfg_mode_best;

        hits_result_t result = search("ACGTT"_dna4, this->index, cfg) | seqan3::views::to<std::vector>;
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    {  // Find best match with a match at the end, allowing a insertion.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::insertion{1}} |
                                          search_cfg_mode_best;

        hits_result_t result = search("ACGT"_dna4, this->index, cfg) | seqan3::views::to<std::vector>;
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    {  // Find best match with a deletion.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::deletion{1}} |
                                          search_cfg_mode_best;

        hits_result_t result = search("AGT"_dna4, this->index, cfg) | seqan3::views::to<std::vector>;
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    { // Find best match with a match at the end, allowing a deletion.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::deletion{1}} |
                                          search_cfg_mode_best;

        hits_result_t result = search("ACGT"_dna4, this->index, cfg) | seqan3::views::to<std::vector>;
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    {  // Find best match with a substitution at the end.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::substitution{1}} |
                                          search_cfg_mode_best;

        hits_result_t result = search("ACGC"_dna4, this->index, cfg) | seqan3::views::to<std::vector>;
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    {  // Find best match with a match at the end, allowing a substitution.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::substitution{1}} |
                                          search_cfg_mode_best;

        hits_result_t result = search("ACGT"_dna4, this->index, cfg) | seqan3::views::to<std::vector>;
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    {  // Find best match with 2 deletions.
        hits_result_t possible_hits2d{{0, 0}, {0, 4}}; // any of 0, 4 ... 1, 5 are not best hits
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::deletion{2}} |
                                          seqan3::search_cfg::mode{seqan3::search_cfg::best};

        hits_result_t result = search("AGTAGT"_dna4, this->index, cfg) | seqan3::views::to<std::vector>;
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits2d.begin(), possible_hits2d.end(), result[0]) != possible_hits2d.end());
    }
}

TYPED_TEST(search_test, search_strategy_all_best)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    {
        auto search_cfg_mode_all_best = seqan3::search_cfg::mode{seqan3::search_cfg::all_best};
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}} |
                                          search_cfg_mode_all_best;

        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8}})); // 1, 5, 9 are not best hits

        EXPECT_EQ(search("AAAA"_dna4, this->index, cfg) | seqan3::views::to<std::vector>, (hits_result_t{})); // no hit
    }
}

TYPED_TEST(search_test, search_strategy_strata)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}} |
                                          seqan3::search_cfg::mode{seqan3::search_cfg::strata{0}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 4}, {0, 8}}));
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}} |
                                          seqan3::search_cfg::mode{seqan3::search_cfg::strata{1}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{{0, 0}, {0, 1}, {0, 4},
                                                                                  {0, 5}, {0, 8}, {0, 9}}));
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}} |
                                          seqan3::search_cfg::mode{seqan3::search_cfg::strata{1}};
        EXPECT_EQ(search("AAAA"_dna4, this->index, cfg) | seqan3::views::to<std::vector>, (hits_result_t{})); // no hit
    }

    // {
    //     // best hit ACGT with 1 error, i.e. 1+1
    //     seqan3::configuration const cfg = seqan3::search_cfg::max_total_error{1} |
    //                                       seqan3::search_cfg::strategy_strata{1};
    //     EXPECT_EQ(uniquify(search(this->index, "CCGT"_dna4, cfg)), (hits_result_t{{0, 0}, {0, 1}, {0, 2}, {0, 3},
    //                                                                               {0, 4}, {0, 5}, {0, 6}, {0, 7},
    //                                                                               {0, 8}, {0, 9}, {0, 1}, {0, 0}}));
    // }

    // {
    //     // best hit ACGT with 1 error, i.e. 1+1
    //     seqan3::configuration const cfg = seqan3::search_cfg::max_total_error{1} |
    //                                       seqan3::search_cfg::strategy_strata{1};
    //     EXPECT_EQ(uniquify(search(this->index, "CCGT"_dna4, cfg)), (hits_result_t{{0, 0}, {0, 1}, {0, 2}, {0, 3},
    //                                                                               {0, 4}, {0, 5}, {0, 6}, {0, 7},
    //                                                                               {0, 8}, {0, 9}, {0, 1}, {0, 0}}));
    // }
}

TYPED_TEST(search_string_test, error_free_string)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search("at"s, this->index)), (hits_result_t{{0, 14}, {0, 18}}));
        EXPECT_EQ(uniquify(search("Jon"s, this->index)), (hits_result_t{}));
    }
}

TYPED_TEST(search_string_test, error_free_raw)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search("at", this->index)), (hits_result_t{{0, 14}, {0, 18}}));
        EXPECT_EQ(uniquify(search("Jon", this->index)), (hits_result_t{}));
    }
}

TYPED_TEST(search_string_test, multiple_queries_string)
{
    using hits_result_t = typename TestFixture::hits_result_t;
    std::vector<std::string> const queries{"at", "Jon"};

    EXPECT_EQ(uniquify(search(queries, this->index)), (hits_result_t{{0, 14}, {0, 18}})); // 2 and 0 hits
}

TYPED_TEST(search_string_test, multiple_queries_raw)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    EXPECT_EQ(uniquify(search({"at", "Jon"}, this->index)), (hits_result_t{{0, 14}, {0, 18}})); // 2 and 0 hits
}

// TYPED_TEST(search_test, return_iterator_index)
// {
// }
//
// TYPED_TEST(search_test, on_hit)
// {
// }
