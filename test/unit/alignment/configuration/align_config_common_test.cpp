// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/band/static_band.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>

using namespace seqan3;

template <typename T>
class alignment_configuration_test : public ::testing::Test
{};

using test_types = ::testing::Types<align_cfg::aligned_ends<std::remove_const_t<decltype(free_ends_all)>>,
                                    align_cfg::band<static_band>,
                                    align_cfg::gap<gap_scheme<>>,
                                    align_cfg::max_error,
                                    align_cfg::mode<detail::global_alignment_type>,
                                    align_cfg::mode<detail::local_alignment_type>,
                                    align_cfg::parallel,
                                    align_cfg::result<>,
                                    align_cfg::scoring<nucleotide_scoring_scheme<int8_t>>>;

TYPED_TEST_CASE(alignment_configuration_test, test_types);

// TODO: this should go to a typed configuration test that also checks the alignment configuration
TEST(alignment_configuration_test, symmetric_configuration)
{
    for (uint8_t i = 0; i < static_cast<uint8_t>(detail::align_config_id::SIZE); ++i)
    {
        // no element can occur twice in a configuration
        EXPECT_FALSE(detail::compatibility_table<detail::align_config_id>[i][i])
            << "There is a TRUE value on the diagonal of the search configuration matrix.";
        for (uint8_t j = 0; j < i; ++j)
        {
            // symmetric matrix
            EXPECT_EQ(detail::compatibility_table<detail::align_config_id>[i][j],
                      detail::compatibility_table<detail::align_config_id>[j][i])
                << "Search configuration matrix is not symmetric.";
        }
    }
}

TEST(alignment_configuration_test, number_of_configs)
{
    // NOTE(rrahn): You must update this test if you add a new value to align_cfg::id
    EXPECT_EQ(static_cast<uint8_t>(detail::align_config_id::SIZE), 11);
}

TYPED_TEST(alignment_configuration_test, config_element)
{
    EXPECT_TRUE((detail::config_element<TypeParam>));
}

TYPED_TEST(alignment_configuration_test, configuration_exists)
{
    configuration cfg{TypeParam{}};
    EXPECT_TRUE(decltype(cfg)::template exists<TypeParam>());
}

template <typename t, typename cfg_t>
void helper_exists()
{
    EXPECT_TRUE(cfg_t::template exists<t>());
}

template <template <typename ...> typename t, typename cfg_t>
void helper_exists()
{
    EXPECT_TRUE(cfg_t::template exists<t>());
}

TYPED_TEST(alignment_configuration_test, configuration_exists_template)
{
    configuration cfg{TypeParam{}};
    helper_exists<TypeParam, decltype(cfg)>();
}
