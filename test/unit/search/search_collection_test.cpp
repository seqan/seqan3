// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <type_traits>

#include <seqan3/range/views/persist.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include <gtest/gtest.h>

#include "helper.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_phred42;

using namespace std::string_literals;

auto ref_id_and_position = seqan3::views::persist
                         | std::views::transform([] (auto && res)
                         {
                             return std::make_pair(res.reference_id(), res.reference_begin_pos());
                         });
auto query_id = seqan3::views::persist | std::views::transform([] (auto && res) { return res.query_id(); });

template <typename index_t>
class search_test : public ::testing::Test
{
    using collection_result_t = std::pair<typename index_t::size_type, typename index_t::size_type>;

public:
    using hits_result_t = std::vector<collection_result_t>;

    std::vector<seqan3::dna4_vector> text{"ACGTACGTACGT"_dna4, "ACGTACGTACGT"_dna4};
    index_t index{text};
};

template <typename index_t>
class search_string_test : public ::testing::Test
{
    using collection_result_t = std::pair<typename index_t::size_type, typename index_t::size_type>;

public:
    using hits_result_t = std::vector<collection_result_t>;

    hits_result_t expected_hits{{0, 14}, {0, 18}, {1, 17}};
    std::vector<std::string> text{"Garfield the fat cat.", "Yet another text at position 1."};
    index_t index{text};
};

using fm_index_types        = ::testing::Types<seqan3::fm_index<seqan3::dna4, seqan3::text_layout::collection>,
                                               seqan3::bi_fm_index<seqan3::dna4, seqan3::text_layout::collection>>;
using fm_index_string_types = ::testing::Types<seqan3::fm_index<char, seqan3::text_layout::collection>,
                                               seqan3::bi_fm_index<char, seqan3::text_layout::collection>>;

TYPED_TEST_SUITE(search_test, fm_index_types, );
TYPED_TEST_SUITE(search_string_test, fm_index_string_types, );

TYPED_TEST(search_test, error_free)
{
    typename TestFixture::hits_result_t expected_hits{{0, 0}, {0, 4}, {0, 8}, {1, 0}, {1, 4}, {1, 8}}; //{refid, pos}
    typename TestFixture::hits_result_t empty_result{};

    {
        // successful and unsuccesful exact search without cfg
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index) | ref_id_and_position, expected_hits);
        EXPECT_RANGE_EQ(search("ACGG"_dna4, this->index) | ref_id_and_position, empty_result);
    }

    {
        // successful and unsuccesful exact search with empty cfg
        seqan3::configuration const cfg;
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | ref_id_and_position, expected_hits);
        EXPECT_RANGE_EQ(search("ACGG"_dna4, this->index, cfg) | ref_id_and_position, empty_result);
    }

    {
        // successful and unsuccesful exact search using empty max_total_error
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{};
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | ref_id_and_position, expected_hits);
        EXPECT_RANGE_EQ(search("ACGG"_dna4, this->index, cfg) | ref_id_and_position, empty_result);
    }

    {
        // successful and unsuccesful exact search using short version of max_total_error
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}};
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | ref_id_and_position, expected_hits);
        EXPECT_RANGE_EQ(search("ACGG"_dna4, this->index, cfg) | ref_id_and_position, empty_result);
    }

    {
        // successful and unsuccesful exact search using max_total_error
        auto const zero_count = seqan3::search_cfg::error_count{0};
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{zero_count} |
                                          seqan3::search_cfg::max_error_substitution{zero_count} |
                                          seqan3::search_cfg::max_error_insertion{zero_count} |
                                          seqan3::search_cfg::max_error_deletion{zero_count};
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | ref_id_and_position, expected_hits);
        EXPECT_RANGE_EQ(search("ACGG"_dna4, this->index, cfg) | ref_id_and_position, empty_result);
    }

    {
        // successful and unsuccesful exact search using empty max_total_error_rate
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_rate{}};
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | ref_id_and_position, expected_hits);
        EXPECT_RANGE_EQ(search("ACGG"_dna4, this->index, cfg) | ref_id_and_position, empty_result);
    }

    {
        // successful and unsuccesful exact search using short version of max_total_error_rate
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_rate{.0}};
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | ref_id_and_position, expected_hits);
        EXPECT_RANGE_EQ(search("ACGG"_dna4, this->index, cfg) | ref_id_and_position, empty_result);
    }

    {
        // successful and unsuccesful exact search using max_total_error_rate
        auto const zero_rate = seqan3::search_cfg::error_rate{.0};
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{zero_rate} |
                                          seqan3::search_cfg::max_error_substitution{zero_rate} |
                                          seqan3::search_cfg::max_error_insertion{zero_rate} |
                                          seqan3::search_cfg::max_error_deletion{zero_rate};
        EXPECT_RANGE_EQ(search("ACGT"_dna4, this->index, cfg) | ref_id_and_position, expected_hits);
        EXPECT_RANGE_EQ(search("ACGG"_dna4, this->index, cfg) | ref_id_and_position, empty_result);
    }
}

TYPED_TEST(search_test, convertible_query)
{
    std::vector<seqan3::qualified<seqan3::dna4, seqan3::phred42>> query{{'A'_dna4, '!'_phred42},
                                                                        {'C'_dna4, '!'_phred42},
                                                                        {'G'_dna4, '!'_phred42},
                                                                        {'T'_dna4, '!'_phred42}};

    typename TestFixture::hits_result_t expected_hits{{0, 0}, {0, 4}, {0, 8}, {1, 0}, {1, 4}, {1, 8}};
    EXPECT_RANGE_EQ(search(query, this->index) | ref_id_and_position, expected_hits);
}

// https://github.com/seqan/seqan3/issues/1542
TYPED_TEST(search_test, single_element_collection)
{
    std::vector<seqan3::dna4_vector> text{"GGTGTGTCGGCCCTATCCCTTGCG"_dna4};
    TypeParam index{text};

    seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{3}} |
                                      seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{3}} |
                                      seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{0}} |
                                      seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{0}};
    typename TestFixture::hits_result_t expected_hits{{0, 0}};
    EXPECT_RANGE_EQ(search("GGTGTGTCG"_dna4, index, cfg) | ref_id_and_position, expected_hits);
}

TYPED_TEST(search_test, multiple_queries)
{
    std::vector<std::vector<seqan3::dna4>> const queries{{"GG"_dna4, "ACGTACGTACGT"_dna4, "ACGTA"_dna4}};

    typename TestFixture::hits_result_t expected_hits{{0, 0}, {1, 0}, {0, 0}, {0, 4}, {1, 0}, {1, 4}};
    std::vector<size_t> expected_query_ids{1, 1, 2, 2, 2, 2};

    seqan3::configuration const cfg = seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_rate{.0}} |
                                      seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_rate{.0}} |
                                      seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_rate{.0}} |
                                      seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_rate{.0}};
    EXPECT_RANGE_EQ(search(queries, this->index, cfg) | ref_id_and_position, expected_hits);
    EXPECT_RANGE_EQ(search(queries, this->index, cfg) | query_id, expected_query_ids);
}

TYPED_TEST(search_string_test, error_free_string)
{
    typename TestFixture::hits_result_t empty_result{};
    // successful and unsuccesful exact search without cfg
    EXPECT_RANGE_EQ(search("at"s, this->index) | ref_id_and_position, this->expected_hits);
    EXPECT_RANGE_EQ(search("Jon"s, this->index) | ref_id_and_position, empty_result);
}

TYPED_TEST(search_string_test, error_free_raw)
{
    typename TestFixture::hits_result_t empty_result{};
    // successful and unsuccesful exact search without cfg
    EXPECT_RANGE_EQ(search("at", this->index) | ref_id_and_position, this->expected_hits);
    EXPECT_RANGE_EQ(search("Jon", this->index) | ref_id_and_position, empty_result);
}

TYPED_TEST(search_string_test, multiple_queries_string)
{
    std::vector<std::string> const queries{"at", "Jon"};

    EXPECT_RANGE_EQ(search(queries, this->index) | ref_id_and_position, this->expected_hits); // 3 and 0 hits
}

TYPED_TEST(search_string_test, multiple_queries_raw)
{
    EXPECT_RANGE_EQ(search({"at", "Jon"}, this->index) | ref_id_and_position, this->expected_hits); // 3 and 0 hits
}
