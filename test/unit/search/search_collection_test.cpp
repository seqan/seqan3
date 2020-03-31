// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <type_traits>

#include <seqan3/search/all.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>

#include <gtest/gtest.h>

#include "helper.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_phred42;

using namespace std::string_literals;

template <typename index_t>
class search_test : public ::testing::Test
{
public:
    using collection_result_t = std::pair<typename index_t::size_type, typename index_t::size_type>;
    using hits_result_t = std::vector<std::pair<size_t, collection_result_t>>;

    std::vector<seqan3::dna4_vector> text{"ACGTACGTACGT"_dna4, "ACGTACGTACGT"_dna4};
    index_t index{text};
};

template <typename index_t>
class search_string_test : public ::testing::Test
{
public:
    using collection_result_t = std::pair<typename index_t::size_type, typename index_t::size_type>;
    using hits_result_t = std::vector<std::pair<size_t, collection_result_t>>;

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
    using hits_result_t = typename TestFixture::hits_result_t;

    hits_result_t expected_hits{{0, {0, 0}}, {0, {0, 4}}, {0, {0, 8}},
                                {0, {1, 0}}, {0, {1, 4}}, {0, {1, 8}}};
    {
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index)), expected_hits);
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index)), hits_result_t{});
    }

    {
        // successful and unsuccesful exact search with empty cfg
        seqan3::configuration const cfg;
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), expected_hits);
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), hits_result_t{});
    }

    {
        // successful and unsuccesful exact search using empty max_total_error
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), expected_hits);
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), hits_result_t{});
    }

    {
        // successful and unsuccesful exact search using short version of max_total_error
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{0}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), expected_hits);
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), hits_result_t{});
    }

    {
        // successful and unsuccesful exact search using max_total_error
        seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{0},
                                                                        seqan3::search_cfg::substitution{0},
                                                                        seqan3::search_cfg::insertion{0},
                                                                        seqan3::search_cfg::deletion{0}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), expected_hits);
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), hits_result_t{});
    }

    {
        // successful and unsuccesful exact search using empty max_total_error_rate
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), expected_hits);
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), hits_result_t{});
    }

    {
        // successful and unsuccesful exact search using short version of max_total_error_rate
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::total{.0}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), expected_hits);
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), hits_result_t{});
    }

    {
        // successful and unsuccesful exact search using max_total_error_rate
        seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::total{.0},
                                                                             seqan3::search_cfg::substitution{.0},
                                                                             seqan3::search_cfg::insertion{.0},
                                                                             seqan3::search_cfg::deletion{.0}};
        EXPECT_EQ(uniquify(search("ACGT"_dna4, this->index, cfg)), expected_hits);
        EXPECT_EQ(uniquify(search("ACGG"_dna4, this->index, cfg)), hits_result_t{});
    }
}

TYPED_TEST(search_test, convertible_query)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    std::vector<seqan3::qualified<seqan3::dna4, seqan3::phred42>> query{{'A'_dna4, '!'_phred42},
                                                                        {'C'_dna4, '!'_phred42},
                                                                        {'G'_dna4, '!'_phred42},
                                                                        {'T'_dna4, '!'_phred42}};

    hits_result_t expected_hits{{0, {0, 0}}, {0, {0, 4}}, {0, {0, 8}},
                                {0, {1, 0}}, {0, {1, 4}}, {0, {1, 8}}};
    EXPECT_EQ(uniquify(search(query, this->index)), expected_hits);
}

TYPED_TEST(search_test, single_element_collection)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    std::vector<seqan3::dna4_vector> text{"ACGATACG"_dna4};
    TypeParam index{text};

    seqan3::configuration const cfg = seqan3::search_cfg::max_error{seqan3::search_cfg::total{1},
                                                                    seqan3::search_cfg::substitution{1},
                                                                    seqan3::search_cfg::insertion{0},
                                                                    seqan3::search_cfg::deletion{0}};
    EXPECT_EQ(uniquify(search("ACGACACG"_dna4, index, cfg)), (hits_result_t{{0, {0, 0}}}));
}

TYPED_TEST(search_test, multiple_queries)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    std::vector<std::vector<seqan3::dna4>> const queries{{"GG"_dna4, "ACGTACGTACGT"_dna4, "ACGTA"_dna4}};

    hits_result_t expected_hits{{1, {0, 0}}, {1, {1, 0}},
                                {2, {0, 0}}, {2, {0, 4}}, {2, {1, 0}}, {2, {1, 4}}};
    seqan3::configuration const cfg = seqan3::search_cfg::max_error_rate{seqan3::search_cfg::total{.0},
                                                                         seqan3::search_cfg::substitution{.0},
                                                                         seqan3::search_cfg::insertion{.0},
                                                                         seqan3::search_cfg::deletion{.0}};
    EXPECT_EQ(uniquify(search(queries, this->index, cfg)), expected_hits);
}

TYPED_TEST(search_string_test, error_free_string)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    {
        hits_result_t expected_hits{{0, {0, 14}}, {0, {0, 18}}, {0, {1, 17}}};
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search("at"s, this->index)), expected_hits);
        EXPECT_EQ(uniquify(search("Jon"s, this->index)), hits_result_t{});
    }
}

TYPED_TEST(search_string_test, error_free_raw)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    {
        hits_result_t expected_hits{{0, {0, 14}}, {0, {0, 18}}, {0, {1, 17}}};
        // successful and unsuccesful exact search without cfg
        EXPECT_EQ(uniquify(search("at", this->index)), expected_hits);
        EXPECT_EQ(uniquify(search("Jon", this->index)), hits_result_t{});
    }
}

TYPED_TEST(search_string_test, multiple_queries_string)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    std::vector<std::string> const queries{"at", "Jon"};

    hits_result_t expected_hits{{0, {0, 14}}, {0, {0, 18}}, {0, {1, 17}}};
    EXPECT_EQ(uniquify(search(queries, this->index)), expected_hits); // 3 and 0 hits
}

TYPED_TEST(search_string_test, multiple_queries_raw)
{
    using hits_result_t = typename TestFixture::hits_result_t;

    hits_result_t expected_hits{{0, {0, 14}}, {0, {0, 18}}, {0, {1, 17}}};
    EXPECT_EQ(uniquify(search({"at", "Jon"}, this->index)), expected_hits); // 3 and 0 hits
}
