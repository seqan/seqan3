// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/search/configuration/hit.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/on_result.hpp>
#include <seqan3/search/configuration/output.hpp>
#include <seqan3/search/configuration/parallel.hpp>
#include <seqan3/search/configuration/result_type.hpp>
#include <seqan3/search/search_result.hpp>

template <typename T>
class search_configuration_test : public ::testing::Test
{};

using search_result_t = seqan3::search_result<seqan3::detail::empty_type,
                                              seqan3::detail::empty_type,
                                              seqan3::detail::empty_type,
                                              seqan3::detail::empty_type>;

// Needed to test the on_result config.
inline constexpr auto on_result_caller = [](auto &&) {};

using test_types = ::testing::Types<seqan3::search_cfg::max_error_total,
                                    seqan3::search_cfg::max_error_substitution,
                                    seqan3::search_cfg::max_error_insertion,
                                    seqan3::search_cfg::max_error_deletion,
                                    seqan3::search_cfg::hit_single_best,
                                    seqan3::search_cfg::on_result<decltype(on_result_caller)>,
                                    seqan3::search_cfg::output_query_id,
                                    seqan3::search_cfg::output_reference_id,
                                    seqan3::search_cfg::output_reference_begin_position,
                                    seqan3::search_cfg::output_index_cursor,
                                    seqan3::search_cfg::parallel,
                                    seqan3::search_cfg::detail::result_type<search_result_t>>;

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

TYPED_TEST(search_configuration_test, config_element)
{
    EXPECT_TRUE((seqan3::detail::config_element<TypeParam>));
}

TYPED_TEST(search_configuration_test, configuration_exists)
{
    seqan3::configuration cfg{TypeParam{}};
    EXPECT_TRUE(decltype(cfg)::template exists<TypeParam>());
}

TEST(search_configuration_test, max_error_defaults)
{
    // empty config defaults to 0 for every empty error configuration
    EXPECT_EQ((seqan3::search_cfg::max_error_total{}.error),
              (seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}}.error));

    EXPECT_EQ((seqan3::search_cfg::max_error_substitution{}.error),
              (seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{0}}.error));

    EXPECT_EQ((seqan3::search_cfg::max_error_insertion{}.error),
              (seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{0}}.error));

    EXPECT_EQ((seqan3::search_cfg::max_error_deletion{}.error),
              (seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{0}}.error));

    // empty config defaults to 0 for every error_count
    EXPECT_EQ((seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{}}.error),
              (seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_count{0}}.error));

    EXPECT_EQ((seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{}}.error),
              (seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_count{0}}.error));

    EXPECT_EQ((seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{}}.error),
              (seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_count{0}}.error));

    EXPECT_EQ((seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{}}.error),
              (seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_count{0}}.error));

    // empty config defaults to 0 for every error_rate
    EXPECT_EQ((seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_rate{}}.error),
              (seqan3::search_cfg::max_error_total{seqan3::search_cfg::error_rate{0}}.error));

    EXPECT_EQ((seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_rate{}}.error),
              (seqan3::search_cfg::max_error_substitution{seqan3::search_cfg::error_rate{0}}.error));

    EXPECT_EQ((seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_rate{}}.error),
              (seqan3::search_cfg::max_error_insertion{seqan3::search_cfg::error_rate{0}}.error));

    EXPECT_EQ((seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_rate{}}.error),
              (seqan3::search_cfg::max_error_deletion{seqan3::search_cfg::error_rate{0}}.error));
}
