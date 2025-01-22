// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <type_traits>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/search/configuration/hit.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/on_result.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "helper.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_phred42;

using namespace std::string_literals;

auto position = std::views::all
              | std::views::transform(
                    [](auto && res)
                    {
                        return res.reference_begin_position();
                    });
auto query_id = std::views::all
              | std::views::transform(
                    [](auto && res)
                    {
                        return res.query_id();
                    });

template <typename index_t>
class search_test : public ::testing::Test
{
public:
    seqan3::dna4_vector text{"ACGTACGTACGT"_dna4};
    index_t index{text};
};

template <typename index_t>
class search_string_test : public ::testing::Test
{
public:
    std::string text{"Garfield the fat cat."};
    index_t index{text};
};

using fm_index_types = ::testing::Types<seqan3::fm_index<seqan3::dna4, seqan3::text_layout::single>,
                                        seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::single>>;
using fm_index_string_types = ::testing::Types<seqan3::fm_index<char, seqan3::text_layout::single>,
                                               seqan3::bi_fm_index<char, seqan3::text_layout::single>>;

TYPED_TEST_SUITE(search_test, fm_index_types, );
TYPED_TEST_SUITE(search_string_test, fm_index_string_types, );

TYPED_TEST(search_test, error_free)
{
    {
        // successful and unsuccesful exact search without cfg
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index) | position, (std::vector{0, 4, 8}));
        EXPECT_RANGE_EQ(search("ACGG"_dna4, this->index) | position, (std::vector<int>{}));
    }

    {
        // successful and unsuccesful exact search with empty cfg
        seqan3::configuration const cfg;
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position, (std::vector{0, 4, 8}));
        EXPECT_RANGE_EQ(search("ACGG"_dna4, this->index, cfg) | position, (std::vector<int>{}));
    }

    {
        // successful and unsuccesful exact search using empty error_count
        // default max_error_total{error_count{}} sets all error to 0
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{}};
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position, (std::vector{0, 4, 8}));
        EXPECT_RANGE_EQ(search("ACGG"_dna4, this->index, cfg) | position, (std::vector<int>{}));
    }

    {
        // successful and unsuccesful exact search using empty error_rate
        // default max_error_total{error_rate{}} sets all error to 0.0
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_rate{}};
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position, (std::vector{0, 4, 8}));
        EXPECT_RANGE_EQ(search("ACGG"_dna4, this->index, cfg) | position, (std::vector<int>{}));
    }
}

TYPED_TEST(search_test, on_result_invocation)
{
    std::vector<int> actual_positions{};
    seqan3::configuration const cfg =
        seqan3::search_cfg::on_result{[&](auto && search_result) -> void
                                      {
                                          actual_positions.push_back(search_result.reference_begin_position());
                                      }};

    search("ACGT"_dna4, this->index, cfg);
    EXPECT_RANGE_EQ(actual_positions, (std::vector{0, 4, 8}));

    actual_positions.clear();
    search("ACGG"_dna4, this->index, cfg);

    EXPECT_TRUE(actual_positions.empty());
}

TYPED_TEST(search_test, convertible_query)
{
    std::vector<seqan3::qualified<seqan3::dna4, seqan3::phred42>> query{{'A'_dna4, '!'_phred42},
                                                                        {'C'_dna4, '!'_phred42},
                                                                        {'G'_dna4, '!'_phred42},
                                                                        {'T'_dna4, '!'_phred42}};

    EXPECT_RANGE_EQ(search(query, this->index) | position, (std::vector{0, 4, 8}));
}

TYPED_TEST(search_test, multiple_queries)
{
    std::vector<std::vector<seqan3::dna4>> const queries{{"GG"_dna4, "ACGTACGTACGT"_dna4, "ACGTA"_dna4}};

    seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_rate{.0}}
                                    | seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_rate{.0}}
                                    | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_rate{.0}}
                                    | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_rate{.0}};

    EXPECT_RANGE_EQ(search(queries, this->index, cfg) | query_id, (std::vector{1, 2, 2}));
    EXPECT_RANGE_EQ(search(queries, this->index, cfg) | position, (std::vector{0, 0, 4}));
}

TYPED_TEST(search_test, parallel_queries)
{
    constexpr size_t num_queries{100u};
    std::vector<std::vector<seqan3::dna4>> const queries{num_queries, {"ACGTACGTACGT"_dna4}};

    seqan3::configuration const cfg =
        seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_rate{.0}}
        | seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_rate{.0}}
        | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_rate{.0}}
        | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_rate{.0}}
        | seqan3::search_cfg::parallel{std::min<uint32_t>(2, std::thread::hardware_concurrency())};

    EXPECT_RANGE_EQ(search(queries, this->index, cfg) | query_id, std::views::iota(0u, num_queries));
    EXPECT_RANGE_EQ(search(queries, this->index, cfg) | position, std::vector(num_queries, 0));
}

TYPED_TEST(search_test, invalid_error_configuration)
{
    seqan3::configuration const cfg1 = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_rate{-0.5}};
    EXPECT_THROW(search("A"_dna4, this->index, cfg1).begin(), std::invalid_argument);

    seqan3::configuration const cfg2 = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_rate{2.0}};
    EXPECT_THROW(search("A"_dna4, this->index, cfg2).begin(), std::invalid_argument);

    seqan3::configuration const cfg3 = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}}
                                     | seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{1}};
    EXPECT_THROW(search("A"_dna4, this->index, cfg3).begin(), std::invalid_argument);

    seqan3::configuration const cfg4 = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}}
                                     | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{1}};
    EXPECT_THROW(search("A"_dna4, this->index, cfg4).begin(), std::invalid_argument);

    seqan3::configuration const cfg5 = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}}
                                     | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{1}};
    EXPECT_THROW(search("A"_dna4, this->index, cfg5).begin(), std::invalid_argument);
}

TYPED_TEST(search_test, invalid_dynamic_hit_configuration)
{
    seqan3::configuration const cfg =
        seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}} | seqan3::search_cfg::hit{};
    EXPECT_THROW(search("A"_dna4, this->index, cfg), std::invalid_argument);
}

TYPED_TEST(search_test, error_substitution)
{
    {
        seqan3::configuration const cfg =
            seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_rate{.25}};

        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position, (std::vector{0, 4, 8})); // exact match
        EXPECT_RANGE_EQ(search("CGG"_dna4, this->index, cfg) | position, (std::vector<int>{})); // not enough mismatches
        EXPECT_RANGE_EQ(search("CGTC"_dna4, this->index, cfg) | position, (std::vector{1, 5})); // 1 mismatch
        EXPECT_RANGE_EQ(search("ACGGACG"_dna4, this->index, cfg) | position, (std::vector{0, 4}));  // 1 mismatch
        EXPECT_RANGE_EQ(search("ACGGACGG"_dna4, this->index, cfg) | position, (std::vector{0, 4})); // 2 mismatches
    }

    {
        seqan3::configuration const cfg =
            seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{1}};

        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position, (std::vector{0, 4, 8})); // exact match
        EXPECT_RANGE_EQ(search("CGTTT"_dna4, this->index, cfg) | position,
                        (std::vector<int>{})); // not enough mismatches
        EXPECT_RANGE_EQ(search("CGG"_dna4, this->index, cfg) | position, (std::vector{1, 5, 9}));  // 1 mismatch
        EXPECT_RANGE_EQ(search("ACGGACG"_dna4, this->index, cfg) | position, (std::vector{0, 4})); // 1 mismatch
        EXPECT_RANGE_EQ(search("CGTCCGTA"_dna4, this->index, cfg) | position, (std::vector{1}));   // 1 mismatch
    }
}

TYPED_TEST(search_test, error_configuration_types)
{
    seqan3::configuration const cfg = seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{1}};
    EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position, (std::vector{0, 4, 8}));
}

TYPED_TEST(search_test, error_insertion)
{
    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_rate{.25}};

        // exact match and insertion at the beginning of the query
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position, (std::vector{0, 1, 4, 5, 8, 9}));
        // 1 insertion
        EXPECT_RANGE_EQ(search("CCGT"_dna4, this->index, cfg) | position, (std::vector{1, 5, 9}));
        // 2 insertions
        EXPECT_RANGE_EQ(search("ACCGGTAC"_dna4, this->index, cfg) | position, (std::vector{0, 4}));
        // 2 insertions necessary, only 1 allowed
        EXPECT_RANGE_EQ(search("ACCGG"_dna4, this->index, cfg) | position, (std::vector<int>{}));
        // deletion necessary, not allowed
        EXPECT_RANGE_EQ(search("ACTACGT"_dna4, this->index, cfg) | position, (std::vector<int>{}));
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{1}};

        // exact match and insertion at the beginning of the query
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position, (std::vector{0, 1, 4, 5, 8, 9}));
        // 1 insertion
        EXPECT_RANGE_EQ(search("CCGT"_dna4, this->index, cfg) | position, (std::vector{1, 5, 9}));
        // 2 insertions necessary, only 1 allowed
        EXPECT_RANGE_EQ(search("ACCGGTAC"_dna4, this->index, cfg) | position, (std::vector<int>{}));
        // deletion necessary, not allowed
        EXPECT_RANGE_EQ(search("ACTACGT"_dna4, this->index, cfg) | position, (std::vector<int>{}));
    }
}

TYPED_TEST(search_test, error_deletion)
{
    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_rate{.25}};

        // exact match, no deletion
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position, (std::vector{0, 4, 8}));
        // not enough max errors
        EXPECT_RANGE_EQ(search("AGT"_dna4, this->index, cfg) | position, (std::vector<int>{}));
        // one deletion (C)
        EXPECT_RANGE_EQ(search("AGTA"_dna4, this->index, cfg) | position, (std::vector{0, 4}));
        // two deletion (C)
        EXPECT_RANGE_EQ(search("AGTAGTAC"_dna4, this->index, cfg) | position, (std::vector{0}));
        // no deletion at beginning. 0 and 4 cannot be reported
        EXPECT_RANGE_EQ(search("CGTACGT"_dna4, this->index, cfg) | position, (std::vector{1, 5}));
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{1}};

        // exact match, no deletion
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position, (std::vector{0, 4, 8}));
        // one deletion (C)
        EXPECT_RANGE_EQ(search("AGTA"_dna4, this->index, cfg) | position, (std::vector{0, 4}));
        // 2 deletions necessary, only 1 allowed
        EXPECT_RANGE_EQ(search("AGTAGTAC"_dna4, this->index, cfg) | position, (std::vector<int>{}));
        // no deletion at beginning. 0 and 4 cannot be reported
        EXPECT_RANGE_EQ(search("CGTACGT"_dna4, this->index, cfg) | position, (std::vector{1, 5}));
    }
}

TYPED_TEST(search_test, error_levenshtein)
{
    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}};
        EXPECT_RANGE_EQ(search("CCGT"_dna4, this->index, cfg) | position, (std::vector{0, 1, 4, 5, 8, 9}));
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{2}};
        EXPECT_RANGE_EQ(search("CCGT"_dna4, this->index, cfg) | position,
                        (std::vector{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10}));
    }
    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{4}};
        EXPECT_RANGE_EQ(search("CCGT"_dna4, this->index, cfg) | position,
                        (std::vector{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}));
    }
}

TYPED_TEST(search_test, error_indel_no_substitution)
{
    {
        // Match one mismatch with 1 insertion and deletion since mismatches are not allowed.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{2}}
                                        | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{2}}
                                        | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{2}};
        EXPECT_RANGE_EQ(search("GTACCTAC"_dna4, this->index, cfg) | position, (std::vector{2}));
    }

    {
        // Enumerate a deletion and match one mismatch with 1 insertion and deletion since mismatches are not allowed.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{3}}
                                        | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{3}}
                                        | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{3}};
        EXPECT_RANGE_EQ(search("GTATCCTAC"_dna4, this->index, cfg) | position, (std::vector{2}));
    }
}

TYPED_TEST(search_test, search_strategy_all)
{
    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}};
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position, (std::vector{0, 1, 4, 5, 8, 9}));
    }

    {
        seqan3::configuration const cfg =
            seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}} | seqan3::search_cfg::hit_all{};
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position, (std::vector{0, 1, 4, 5, 8, 9}));
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::hit{seqan3::search_cfg::hit_all{}};
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position, (std::vector{0, 1, 4, 5, 8, 9}));
    }
}

TYPED_TEST(search_test, search_strategy_best)
{
    std::vector possible_hits{0, 4, 8}; // any of 0, 4, 8 ... 1, 5, 9 are not best hits

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::hit_single_best{};

        std::vector result = search("ACGT"_dna4, this->index, cfg) | position | seqan3::ranges::to<std::vector>();
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());

        EXPECT_RANGE_EQ(search("AAAA"_dna4, this->index, cfg) | position, (std::vector<int>{})); // no hit
    }

    { // Find best match with 1 insertion at the end.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::hit_single_best{};

        std::vector result = search("ACGTT"_dna4, this->index, cfg) | position | seqan3::ranges::to<std::vector>();
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    { // Find best match with a match at the end, allowing a insertion.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::hit_single_best{};

        std::vector result = search("ACGT"_dna4, this->index, cfg) | position | seqan3::ranges::to<std::vector>();
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    { // Find best match with a deletion.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::hit_single_best{};

        std::vector result = search("AGT"_dna4, this->index, cfg) | position | seqan3::ranges::to<std::vector>();
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    { // Find best match with a match at the end, allowing a deletion.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::hit_single_best{};

        std::vector result = search("ACGT"_dna4, this->index, cfg) | position | seqan3::ranges::to<std::vector>();
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    { // Find best match with a substitution at the end.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::hit_single_best{};

        std::vector result = search("ACGC"_dna4, this->index, cfg) | position | seqan3::ranges::to<std::vector>();
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    { // Find best match with a match at the end, allowing a substitution.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::hit_single_best{};

        std::vector result = search("ACGT"_dna4, this->index, cfg) | position | seqan3::ranges::to<std::vector>();
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits.begin(), possible_hits.end(), result[0]) != possible_hits.end());
    }

    {                                      // Find best match with 2 deletions.
        std::vector possible_hits2d{0, 4}; // any of 0, 4 ... 1, 5 are not best hits
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{2}}
                                        | seqan3::search_cfg::hit_single_best{};

        std::vector result = search("AGTAGT"_dna4, this->index, cfg) | position | seqan3::ranges::to<std::vector>();
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits2d.begin(), possible_hits2d.end(), result[0]) != possible_hits2d.end());
    }

    {                                      // Find best match with 2 deletions.
        std::vector possible_hits2d{0, 4}; // any of 0, 4 ... 1, 5 are not best hits
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{2}}
                                        | seqan3::search_cfg::hit{seqan3::search_cfg::hit_single_best{}};

        std::vector result = search("AGTAGT"_dna4, this->index, cfg) | position | seqan3::ranges::to<std::vector>();
        ASSERT_EQ(result.size(), 1u);
        EXPECT_TRUE(std::find(possible_hits2d.begin(), possible_hits2d.end(), result[0]) != possible_hits2d.end());
    }

    { // Find best match with at most 2 of every error type.
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{2}}
                                        | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{2}}
                                        | seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{2}}
                                        | seqan3::search_cfg::hit_single_best{};

        EXPECT_RANGE_EQ(search("TATTA"_dna4, this->index, cfg) | position, (std::vector{3}));
    }
}

TYPED_TEST(search_test, search_strategy_all_best)
{
    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::hit_all_best{};

        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position,
                        (std::vector{0, 4, 8})); // 1, 5, 9 are not best hits

        EXPECT_RANGE_EQ(search("AAAA"_dna4, this->index, cfg) | position, (std::vector<int>{})); // no hit
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::hit{seqan3::search_cfg::hit_all_best{}};

        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position,
                        (std::vector{0, 4, 8})); // 1, 5, 9 are not best hits

        EXPECT_RANGE_EQ(search("AAAA"_dna4, this->index, cfg) | position, (std::vector<int>{})); // no hit
    }
}

TYPED_TEST(search_test, search_strategy_strata)
{
    {
        seqan3::configuration const cfg =
            seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}} | seqan3::search_cfg::hit_strata{0};
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position, (std::vector{0, 4, 8}));
    }

    {
        seqan3::configuration const cfg =
            seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}} | seqan3::search_cfg::hit_strata{1};
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | position, (std::vector{0, 1, 4, 5, 8, 9}));
    }

    {
        seqan3::configuration const cfg =
            seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}} | seqan3::search_cfg::hit_strata{1};
        EXPECT_RANGE_EQ(search("AAAA"_dna4, this->index, cfg) | position, (std::vector<int>{})); // no hit
    }

    {
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}}
                                        | seqan3::search_cfg::hit{seqan3::search_cfg::hit_strata{1}};
        EXPECT_RANGE_EQ(search("AAAA"_dna4, this->index, cfg) | position, (std::vector<int>{})); // no hit
    }

    // {
    //     // best hit ACGT with 1 error, i.e. 1+1
    //     seqan3::configuration const cfg = seqan3::search_cfg::max_total_error{1} |
    //                                       seqan3::search_cfg::hit_strata{1};
    //     EXPECT_RANGE_EQ(search(this->index, "CCGT"_dna4, cfg) | position, (std::vector{0, 1, 2, 3,4, 5, 6, 7,
    //                                                                                    8, 9, 1, 0}));
    // }

    // {
    //     // best hit ACGT with 1 error, i.e. 1+1
    //     seqan3::configuration const cfg = seqan3::search_cfg::max_total_error{1} |
    //                                       seqan3::search_cfg::hit_strata{1};
    //     EXPECT_RANGE_EQ(search(this->index, "CCGT"_dna4, cfg) | position, (std::vector{0, 1, 2, 3, 4, 5, 6, 7,
    //                                                                                    8, 9, 1, 0}));
    // }
}

TYPED_TEST(search_test, parallel_without_parameter)
{
    seqan3::configuration cfg = seqan3::search_cfg::parallel{};

    EXPECT_THROW(search("AAAA"_dna4, this->index, cfg), std::runtime_error);
}

TYPED_TEST(search_test, debug_streaming)
{
    std::ostringstream oss;
    seqan3::debug_stream_type stream{oss};
    stream << search("TAC"_dna4, this->index);
    EXPECT_EQ(oss.str(),
              "[<query_id:0, reference_id:0, reference_pos:3>"
              ",<query_id:0, reference_id:0, reference_pos:7>]");
}

// https://github.com/seqan/seqan3/issues/2115
TYPED_TEST(search_test, issue_2115)
{
    std::vector<seqan3::dna4> const genome{"ACAG"_dna4};
    TypeParam const index{genome};

    // One substitution error, report all best hits
    seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{1}}
                                    | seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{1}}
                                    | seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{0}}
                                    | seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{0}}
                                    | seqan3::search_cfg::hit_all_best{};

    std::vector<seqan3::dna4> const dna4_query{"ACGG"_dna4};
    std::vector<seqan3::qualified<seqan3::dna4, seqan3::phred42>> const dna4q_query{{'A'_dna4, '!'_phred42},
                                                                                    {'C'_dna4, '0'_phred42},
                                                                                    {'G'_dna4, '?'_phred42},
                                                                                    {'G'_dna4, 'J'_phred42}};

    // Quality should not alter search results.
    EXPECT_RANGE_EQ(seqan3::search(dna4q_query, index, cfg), seqan3::search(dna4_query, index, cfg));
}

TYPED_TEST(search_string_test, error_free_string)
{
    // successful and unsuccesful exact search without cfg
    EXPECT_RANGE_EQ(search("at"s, this->index) | position, (std::vector{14, 18}));
    EXPECT_RANGE_EQ(search("Jon"s, this->index) | position, (std::vector<int>{}));
}

TYPED_TEST(search_string_test, error_free_raw)
{
    // successful and unsuccesful exact search without cfg
    EXPECT_RANGE_EQ(search("at", this->index) | position, (std::vector{14, 18}));
    EXPECT_RANGE_EQ(search("Jon", this->index) | position, (std::vector<int>{}));
}

TYPED_TEST(search_string_test, multiple_queries_string)
{
    std::vector<std::string> const queries{"at", "Jon"};

    EXPECT_RANGE_EQ(search(queries, this->index) | position, (std::vector{14, 18})); // 2 and 0 hits
}

TYPED_TEST(search_string_test, multiple_queries_raw)
{
    EXPECT_RANGE_EQ(search({"at", "Jon"}, this->index) | position, (std::vector{14, 18})); // 2 and 0 hits
}

// TYPED_TEST(search_test, return_iterator_index)
// {
// }
//
// TYPED_TEST(search_test, on_hit)
// {
// }
