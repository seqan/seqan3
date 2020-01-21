// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <type_traits>

#include <seqan3/search/algorithm/all.hpp>

#include <gtest/gtest.h>

using namespace seqan3;

template <typename T>
class search_configuration_test : public ::testing::Test
{};

using test_types = ::testing::Types<search_cfg::max_error_rate<>,
                                    search_cfg::max_error<>,
                                    search_cfg::mode<detail::search_mode_best>,
                                    search_cfg::output<detail::search_output_text_position>,
                                    search_cfg::parallel>;

TYPED_TEST_SUITE(search_configuration_test, test_types, );

// TODO: this should go to a typed configuration test that also checks the alignment configuration
TEST(search_configuration_test, symmetric_configuration)
{
    for (uint8_t i = 0; i < static_cast<uint8_t>(detail::search_config_id::SIZE); ++i)
    {
        // no element can occur twice in a configuration
        EXPECT_FALSE(detail::compatibility_table<detail::search_config_id>[i][i])
            << "There is a TRUE value on the diagonal of the search configuration matrix.";
        for (uint8_t j = 0; j < i; ++j)
        {
            // symmetric matrix
            EXPECT_EQ(detail::compatibility_table<detail::search_config_id>[i][j],
                      detail::compatibility_table<detail::search_config_id>[j][i])
                << "Search configuration matrix is not symmetric.";
        }
    }
}

TYPED_TEST(search_configuration_test, config_element)
{
    EXPECT_TRUE((detail::config_element<TypeParam>));
}

TYPED_TEST(search_configuration_test, configuration_exists)
{
    configuration cfg{TypeParam{}};
    EXPECT_TRUE(decltype(cfg)::template exists<TypeParam>());
}

// TEST(search_configuration_test, illegal_runtime_configurations)
// {
//     std::vector<dna4> text{"ACGT"_dna4}, query{"ACG"_dna4};
//     fm_index<std::vector<dna4>> fm{text};
//
//     // max_error* without error_type
//     search(fm, query, max_total_error(0) | error_type(error_type_enum::none));
//     EXPECT_DEATH(search(fm, query, detail::configuration {max_error(1)}), "");
//     EXPECT_DEATH(search(fm, query, detail::configuration {max_error_rate(.1)}), "");
//     EXPECT_DEATH(search(fm, query, max_error(1) | error_type(error_type_enum::none)), "");
//
//     // error_type without max_error*
//     search(fm, query, detail::configuration {error_type(error_type_enum::none)});
//     search(fm, query, error_type(error_type_enum::none) | max_error_rate(.0));
//     EXPECT_DEATH(search(fm, query, detail::configuration {error_type(error_type_enum::substitution)}), "");
//     EXPECT_DEATH(search(fm, query, error_type(error_type_enum::substitution) | max_error(0)), "");
//     EXPECT_DEATH(search(fm, query, error_type(error_type_enum::substitution) | max_error_rate(.0)), "");
//
//     // error_type with max_error*
//     search(fm, query, error_type(error_type_enum::substitution) | max_error(1));
//     search(fm, query, error_type(error_type_enum::substitution) | max_error_rate(.1));
// }
