// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include <seqan3/search/configuration/hit.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/output.hpp>
#include <seqan3/search/configuration/parallel.hpp>
#include <seqan3/search/configuration/result_type.hpp>
#include <seqan3/search/search_result.hpp>

#include <gtest/gtest.h>

template <typename T>
class search_configuration_test : public ::testing::Test
{};

using search_result_t = seqan3::search_result<seqan3::detail::empty_type,
                                              seqan3::detail::empty_type,
                                              seqan3::detail::empty_type,
                                              seqan3::detail::empty_type>;

using test_types = ::testing::Types<seqan3::search_cfg::max_error_total,
                                    seqan3::search_cfg::max_error_substitution,
                                    seqan3::search_cfg::max_error_insertion,
                                    seqan3::search_cfg::max_error_deletion,
                                    seqan3::detail::hit_single_best_tag,
                                    seqan3::detail::output_query_id_tag,
                                    seqan3::detail::output_reference_id_tag,
                                    seqan3::detail::output_reference_begin_position_tag,
                                    seqan3::detail::output_index_cursor_tag,
                                    seqan3::search_cfg::parallel,
                                    seqan3::search_cfg::detail::result_type_tag<search_result_t>>;

TYPED_TEST_SUITE(search_configuration_test, test_types, );

// TODO: this should go to a typed configuration test that also checks the alignment configuration
TEST(search_configuration_test, symmetric_configuration)
{
    for (uint8_t i = 0; i < static_cast<uint8_t>(seqan3::detail::search_config_id::SIZE); ++i)
    {
        // no element can occur twice in a configuration
        EXPECT_FALSE(seqan3::detail::compatibility_table<seqan3::detail::search_config_id>[i][i])
            << "There is a TRUE value on the diagonal of the search configuration matrix.";
        for (uint8_t j = 0; j < i; ++j)
        {
            // symmetric matrix
            EXPECT_EQ(seqan3::detail::compatibility_table<seqan3::detail::search_config_id>[i][j],
                      seqan3::detail::compatibility_table<seqan3::detail::search_config_id>[j][i])
                << "Search configuration matrix is not symmetric.";
        }
    }
}

TYPED_TEST(search_configuration_test, config_element_specialisation)
{
    EXPECT_TRUE((seqan3::detail::config_element_specialisation<TypeParam>));
}

TYPED_TEST(search_configuration_test, configuration_exists)
{
    seqan3::configuration cfg{TypeParam{}};
    EXPECT_TRUE(decltype(cfg)::template exists<TypeParam>());
}

TEST(search_configuration_test, max_error_defaults)
{
    // empty config defaults to 0 for every empty error configuration
    EXPECT_EQ((seqan3::search_cfg::max_error_total{}.value),
              (seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}}.value));

    EXPECT_EQ((seqan3::search_cfg::max_error_substitution{}.value),
              (seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{0}}.value));

    EXPECT_EQ((seqan3::search_cfg::max_error_insertion{}.value),
              (seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{0}}.value));

    EXPECT_EQ((seqan3::search_cfg::max_error_deletion{}.value),
              (seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{0}}.value));

    // empty config defaults to 0 for every error_count
    EXPECT_EQ((seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{}}.value),
              (seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}}.value));

    EXPECT_EQ((seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{}}.value),
              (seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{0}}.value));

    EXPECT_EQ((seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{}}.value),
              (seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{0}}.value));

    EXPECT_EQ((seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{}}.value),
              (seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{0}}.value));

    // empty config defaults to 0 for every error_rate
    EXPECT_EQ((seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_rate{}}.value),
              (seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_rate{0}}.value));

    EXPECT_EQ((seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_rate{}}.value),
              (seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_rate{0}}.value));

    EXPECT_EQ((seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_rate{}}.value),
              (seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_rate{0}}.value));

    EXPECT_EQ((seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_rate{}}.value),
              (seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_rate{0}}.value));
}
