// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <type_traits>

#include <seqan3/search/algorithm/all.hpp>

#include <gtest/gtest.h>

#include "helper.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_phred42;

using namespace std::string_literals;

template <typename T>
class search_test : public ::testing::Test
{
public:
    seqan3::dna4_vector text{"ACGTACGTACGT"_dna4};
    T index{text};
};

template <typename T>
class search_string_test : public ::testing::Test
{
public:
    std::string text{"Garfield the fat cat."};
    T index{text};
};

using fm_index_types        = ::testing::Types<seqan3::fm_index<seqan3::dna4, seqan3::text_layout::single>,
                                               seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::single>>;
using fm_index_string_types = ::testing::Types<seqan3::fm_index<char, seqan3::text_layout::single>,
                                               seqan3::bi_fm_index<char, seqan3::text_layout::single>>;

TYPED_TEST_SUITE(search_test, fm_index_types, );
TYPED_TEST_SUITE(search_string_test, fm_index_string_types, );

TYPED_TEST(search_test, error_free)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4, this->index)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(seqan3::uniquify(search("ACGG"_dna4, this->index)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search with empty cfg
        seqan3::configuration const cfg;
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(seqan3::uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using empty max_total_error
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{};
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(seqan3::uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using short version of max_total_error
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}};
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(seqan3::uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using max_total_error
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{0},
                                                                        seqan3::search_cfg::substitution{0},
                                                                        seqan3::search_cfg::insertion{0},
                                                                        seqan3::search_cfg::deletion{0}};
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(seqan3::uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using empty max_total_error_rate
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{};
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(seqan3::uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using short version of max_total_error_rate
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::total{.0}};
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(seqan3::uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }

    {
        // successful and unsuccesful exact search using max_total_error_rate
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::total{.0},
                                                                             seqan3::search_cfg::substitution{.0},
                                                                             seqan3::search_cfg::insertion{.0},
                                                                             seqan3::search_cfg::deletion{.0}};
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
        EXPECT_EQ(seqan3::uniquify(search("ACGG"_dna4, this->index, cfg)), (hits_result_t{}));
    }
}

TYPED_TEST(search_test, convertible_query)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;
    std::vector<seqan3::qualified<seqan3::dna4, seqan3::phred42>> query{{'A'_dna4, '!'_phred42},
                                                                        {'C'_dna4, '!'_phred42},
                                                                        {'G'_dna4, '!'_phred42},
                                                                        {'T'_dna4, '!'_phred42}};

    EXPECT_EQ(seqan3::uniquify(search(query, this->index)), (hits_result_t{0, 4, 8}));
}

TYPED_TEST(search_test, multiple_queries)
{
    using hits_result_t = std::vector<std::vector<typename TypeParam::size_type>>;
    std::vector<seqan3::dna4_vector> const queries{{"GG"_dna4, "ACGTACGTACGT"_dna4, "ACGTA"_dna4}};

    seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::total{.0},
                                                                         seqan3::search_cfg::substitution{.0},
                                                                         seqan3::search_cfg::insertion{.0},
                                                                         seqan3::search_cfg::deletion{.0}};
    EXPECT_EQ(seqan3::uniquify(search(queries, this->index, cfg)), (hits_result_t{{}, {0}, {0, 4}})); // 0, 1 and 2 hits
}

TYPED_TEST(search_test, invalid_error_configuration)
{
    seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{0},
                                                                    seqan3::search_cfg::substitution{1}};
    EXPECT_THROW(search("A"_dna4, this->index, cfg), std::invalid_argument);
}

TYPED_TEST(search_test, error_substitution)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::total{.25},
                                                                             seqan3::search_cfg::substitution{.25}};

        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4,     this->index, cfg)), (hits_result_t{0, 4, 8})); // exact match
        EXPECT_EQ(seqan3::uniquify(search("CGG"_dna4,      this->index, cfg)), (hits_result_t{}));        // not enough
                                                                                                          //  mismatches
        EXPECT_EQ(seqan3::uniquify(search("CGTC"_dna4,     this->index, cfg)), (hits_result_t{1, 5}));    // 1 mismatch
        EXPECT_EQ(seqan3::uniquify(search("ACGGACG"_dna4,  this->index, cfg)), (hits_result_t{0, 4}));    // 1 mismatch
        EXPECT_EQ(seqan3::uniquify(search("ACGGACGG"_dna4, this->index, cfg)), (hits_result_t{0, 4}));    // 2 mismatches
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::total{.25},
                                                                             seqan3::search_cfg::substitution{.25},
                                                                             seqan3::search_cfg::insertion{.0},
                                                                             seqan3::search_cfg::deletion{.0}};

        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4,     this->index, cfg)), (hits_result_t{0, 4, 8})); // exact match
        EXPECT_EQ(seqan3::uniquify(search("CGG"_dna4,      this->index, cfg)), (hits_result_t{}));        // not enough
                                                                                                          //  mismatches
        EXPECT_EQ(seqan3::uniquify(search("CGTC"_dna4,     this->index, cfg)), (hits_result_t{1, 5}));    // 1 mismatch
        EXPECT_EQ(seqan3::uniquify(search("ACGGACG"_dna4,  this->index, cfg)), (hits_result_t{0, 4}));    // 1 mismatch
        EXPECT_EQ(seqan3::uniquify(search("ACGGACGG"_dna4, this->index, cfg)), (hits_result_t{0, 4}));    // 2 mismatches
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1},
                                                                        seqan3::search_cfg::substitution{1}};

        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4,     this->index, cfg)), (hits_result_t{0, 4, 8})); // exact match
        EXPECT_EQ(seqan3::uniquify(search("CGTTT"_dna4,    this->index, cfg)), (hits_result_t{}));        // not enough
                                                                                                          //  mismatches
        EXPECT_EQ(seqan3::uniquify(search("CGG"_dna4,      this->index, cfg)), (hits_result_t{1, 5, 9})); // 1 mismatch
        EXPECT_EQ(seqan3::uniquify(search("ACGGACG"_dna4,  this->index, cfg)), (hits_result_t{0, 4}));    // 1 mismatch
        EXPECT_EQ(seqan3::uniquify(search("CGTCCGTA"_dna4, this->index, cfg)), (hits_result_t{1}));       // 1 mismatch
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1},
                                                                        seqan3::search_cfg::substitution{1},
                                                                        seqan3::search_cfg::insertion{0},
                                                                        seqan3::search_cfg::deletion{0}};

        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4,     this->index, cfg)), (hits_result_t{0, 4, 8})); // exact match
        EXPECT_EQ(seqan3::uniquify(search("CGTTT"_dna4,    this->index, cfg)), (hits_result_t{}));        // not enough
                                                                                                          //  mismatches
        EXPECT_EQ(seqan3::uniquify(search("CGG"_dna4,      this->index, cfg)), (hits_result_t{1, 5, 9})); // 1 mismatch
        EXPECT_EQ(seqan3::uniquify(search("ACGGACG"_dna4,  this->index, cfg)), (hits_result_t{0, 4}));    // 1 mismatch
        EXPECT_EQ(seqan3::uniquify(search("CGTCCGTA"_dna4, this->index, cfg)), (hits_result_t{1}));       // 1 mismatch
    }
}

TYPED_TEST(search_test, error_configuration_types)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        uint8_t s = 1, t = 1;
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{t},
                                                                        seqan3::search_cfg::substitution{s}};
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::substitution{1}};
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
    }
}

TYPED_TEST(search_test, error_insertion)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::total{.25},
                                                                             seqan3::search_cfg::insertion{.25}};

        // exact match and insertion at the beginning of the query
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4,     this->index, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
        // 1 insertion
        EXPECT_EQ(seqan3::uniquify(search("CCGT"_dna4,     this->index, cfg)), (hits_result_t{1, 5, 9}));
        // 2 insertions
        EXPECT_EQ(seqan3::uniquify(search("ACCGGTAC"_dna4, this->index, cfg)), (hits_result_t{0, 4}));
        // 2 insertions necessary, only 1 allowed
        EXPECT_EQ(seqan3::uniquify(search("ACCGG"_dna4,    this->index, cfg)), (hits_result_t{}));
        // deletion necessary, not allowed
        EXPECT_EQ(seqan3::uniquify(search("ACTACGT"_dna4,  this->index, cfg)), (hits_result_t{}));
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1},
                                                                        seqan3::search_cfg::insertion{1}};

        // exact match and insertion at the beginning of the query
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4,     this->index, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
        // 1 insertion
        EXPECT_EQ(seqan3::uniquify(search("CCGT"_dna4,     this->index, cfg)), (hits_result_t{1, 5, 9}));
        // 2 insertions necessary, only 1 allowed
        EXPECT_EQ(seqan3::uniquify(search("ACCGGTAC"_dna4, this->index, cfg)), (hits_result_t{}));
        // deletion necessary, not allowed
        EXPECT_EQ(seqan3::uniquify(search("ACTACGT"_dna4,  this->index, cfg)), (hits_result_t{}));
    }
}

TYPED_TEST(search_test, error_deletion)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::total{.25},
                                                                             seqan3::search_cfg::deletion{.25}};

        // exact match, no deletion
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4,     this->index, cfg)), (hits_result_t{0, 4, 8}));
        // not enough max errors
        EXPECT_EQ(seqan3::uniquify(search("AGT"_dna4,      this->index, cfg)), (hits_result_t{}));
        // one deletion (C)
        EXPECT_EQ(seqan3::uniquify(search("AGTA"_dna4,     this->index, cfg)), (hits_result_t{0, 4}));
        // two deletion (C)
        EXPECT_EQ(seqan3::uniquify(search("AGTAGTAC"_dna4, this->index, cfg)), (hits_result_t{0}));
        // no deletion at beginning. 0 and 4 cannot be reported
        EXPECT_EQ(seqan3::uniquify(search("CGTACGT"_dna4,  this->index, cfg)), (hits_result_t{1, 5}));
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1},
                                                                        seqan3::search_cfg::deletion{1}};

        // exact match, no deletion
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4,     this->index, cfg)), (hits_result_t{0, 4, 8}));
        // one deletion (C)
        EXPECT_EQ(seqan3::uniquify(search("AGTA"_dna4,     this->index, cfg)), (hits_result_t{0, 4}));
        // 2 deletions necessary, only 1 allowed
        EXPECT_EQ(seqan3::uniquify(search("AGTAGTAC"_dna4, this->index, cfg)), (hits_result_t{}));
        // no deletion at beginning. 0 and 4 cannot be reported
        EXPECT_EQ(seqan3::uniquify(search("CGTACGT"_dna4,  this->index, cfg)), (hits_result_t{1, 5}));
    }
}

TYPED_TEST(search_test, error_levenshtein)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}};
        EXPECT_EQ(seqan3::uniquify(search("CCGT"_dna4, this->index, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{2}};
        EXPECT_EQ(seqan3::uniquify(search("CCGT"_dna4, this->index, cfg)),
                  (hits_result_t{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
    }
}

TYPED_TEST(search_test, error_indel_no_substitution)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        // Match one mismatch with 1 insertion and deletion since mismatches are not allowed.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{2},
                                                                        seqan3::search_cfg::deletion{2},
                                                                        seqan3::search_cfg::insertion{2}};
        EXPECT_EQ(seqan3::uniquify(search("GTACCTAC"_dna4, this->index, cfg)), (hits_result_t{2}));
    }

    {
        // Enumerate a deletion and match one mismatch with 1 insertion and deletion since mismatches are not allowed.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{3},
                                                                        seqan3::search_cfg::deletion{3},
                                                                        seqan3::search_cfg::insertion{3}};
        EXPECT_EQ(seqan3::uniquify(search("GTATCCTAC"_dna4, this->index, cfg)), (hits_result_t{2}));
    }
}

TYPED_TEST(search_test, search_strategy_all)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}};
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }

    {
        auto search_cfg_mode_all = seqan3::search_cfg::mode{seqan3::search_cfg::all};

        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}}
                                                                        | search_cfg_mode_all;
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }
}

TYPED_TEST(search_test, search_strategy_best)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;
    auto search_cfg_mode_best = seqan3::search_cfg::mode{seqan3::search_cfg::best};

    hits_result_t possible_hits{0, 4, 8}; // any of 0, 4, 8 ... 1, 5, 9 are not best hits
    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}}
                                                                        | search_cfg_mode_best;

        hits_result_t result = search("ACGT"_dna4, this->index, cfg);
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());

        EXPECT_EQ(search("AAAA"_dna4, this->index, cfg), (hits_result_t{})); // no hit
    }

    { // Find best match with 1 insertion at the end.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1},
                                                                        seqan3::search_cfg::insertion{1}}
                                                                        | search_cfg_mode_best;

        hits_result_t result = search("ACGTT"_dna4, this->index, cfg);
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    {  // Find best match with a match at the end, allowing a insertion.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1},
                                                                        seqan3::search_cfg::insertion{1}}
                                                                        | search_cfg_mode_best;

        hits_result_t result = search("ACGT"_dna4, this->index, cfg);
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    {  // Find best match with a deletion.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1},
                                                                        seqan3::search_cfg::deletion{1}}
                                                                        | search_cfg_mode_best;

        hits_result_t result = search("AGT"_dna4, this->index, cfg);
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    { // Find best match with a match at the end, allowing a deletion.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1},
                                                                        seqan3::search_cfg::deletion{1}}
                                                                        | search_cfg_mode_best;

        hits_result_t result = search("ACGT"_dna4, this->index, cfg);
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    {  // Find best match with a substitution at the end.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1},
                                                                        seqan3::search_cfg::substitution{1}}
                                                                        | search_cfg_mode_best;

        hits_result_t result = search("ACGC"_dna4, this->index, cfg);
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    {  // Find best match with a match at the end, allowing a substitution.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1},
                                                                        seqan3::search_cfg::substitution{1}}
                                                                        | search_cfg_mode_best;

        hits_result_t result = search("ACGT"_dna4, this->index, cfg);
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    {  // Find best match with 2 deletions.
        hits_result_t possible_hits2d{0, 4}; // any of 0, 4 ... 1, 5 are not best hits
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{2},
                                                                        seqan3::search_cfg::deletion{2}}
                                                                        | search_cfg_mode_best;

        hits_result_t result = search("AGTAGT"_dna4, this->index, cfg);
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits2d.begin(), possible_hits2d.end(), result[0]) != possible_hits2d.end());
    }
}

TYPED_TEST(search_test, search_strategy_all_best)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        auto search_cfg_mode_all_best = seqan3::search_cfg::mode{seqan3::search_cfg::all_best};
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}}
                                                                        | search_cfg_mode_all_best;

        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4, this->index, cfg)),
                                   (hits_result_t{0, 4, 8})); // 1, 5, 9 are not best hits

        EXPECT_EQ(search("AAAA"_dna4, this->index, cfg), (hits_result_t{})); // no hit
    }
}

TYPED_TEST(search_test, search_strategy_strata)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        auto search_cfg_mode_strata_0 = seqan3::search_cfg::mode{seqan3::search_cfg::strata{0}};
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}}
                                                                        | search_cfg_mode_strata_0;
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 4, 8}));
    }

    {
        auto search_cfg_mode_strata_1 = seqan3::search_cfg::mode{seqan3::search_cfg::strata{1}};
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}}
                                                                        | search_cfg_mode_strata_1;
        EXPECT_EQ(seqan3::uniquify(search("ACGT"_dna4, this->index, cfg)), (hits_result_t{0, 1, 4, 5, 8, 9}));
    }

    {
        auto search_cfg_mode_strata_1 = seqan3::search_cfg::mode{seqan3::search_cfg::strata{1}};
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1}}
                                                                        | search_cfg_mode_strata_1;
        EXPECT_EQ(search("AAAA"_dna4, this->index, cfg), (hits_result_t{})); // no hit
    }

    // {
    //     // best hit ACGT with 1 error, i.e. 1+1
    //     seqan3::configuration const cfg = max_total_error{1} | strategy_strata{1};
    //     EXPECT_EQ(seqan3::uniquify(search(this->index, "CCGT"_dna4, cfg)),
    //                                (hits_result_t{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
    // }

    // {
    //     // best hit ACGT with 1 error, i.e. 1+1
    //     seqan3::configuration const cfg = max_total_error{1} | strategy_strata{1};
    //     EXPECT_EQ(seqan3::uniquify(search(this->index, "CCGT"_dna4, cfg)),
    //                                (hits_result_t{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
    // }
}

TYPED_TEST(search_string_test, error_free_string)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(seqan3::uniquify(search("at"s, this->index)), (hits_result_t{14, 18}));
        EXPECT_EQ(seqan3::uniquify(search("Jon"s, this->index)), (hits_result_t{}));
    }
}

TYPED_TEST(search_string_test, error_free_raw)
{
    using hits_result_t = std::vector<typename TypeParam::size_type>;

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(seqan3::uniquify(search("at", this->index)), (hits_result_t{14, 18}));
        EXPECT_EQ(seqan3::uniquify(search("Jon", this->index)), (hits_result_t{}));
    }
}

TYPED_TEST(search_string_test, multiple_queries_string)
{
    using hits_result_t = std::vector<std::vector<typename TypeParam::size_type>>;
    std::vector<std::string> const queries{"at", "Jon"};

    EXPECT_EQ(seqan3::uniquify(search(queries, this->index)), (hits_result_t{{14, 18}, {}})); // 2 and 0 hits
}

TYPED_TEST(search_string_test, multiple_queries_raw)
{
    using hits_result_t = std::vector<std::vector<typename TypeParam::size_type>>;

    EXPECT_EQ(seqan3::uniquify(search({"at", "Jon"}, this->index)), (hits_result_t{{14, 18}, {}})); // 2 and 0 hits
}

// TYPED_TEST(search_test, return_iterator_index)
// {
// }
//
// TYPED_TEST(search_test, on_hit)
// {
// }
