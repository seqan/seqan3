// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/band/static_band.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>

template <typename T>
class alignment_configuration_test : public ::testing::Test
{};

using alignment_result_t = seqan3::alignment_result<seqan3::detail::alignment_result_value_type<int, int>>;

using test_types = ::testing::Types<seqan3::align_cfg::aligned_ends<std::remove_const_t<decltype(seqan3::free_ends_all)>>,
                                    seqan3::align_cfg::band<seqan3::static_band>,
                                    seqan3::align_cfg::gap<seqan3::gap_scheme<>>,
                                    seqan3::align_cfg::max_error,
                                    seqan3::align_cfg::mode<seqan3::detail::global_alignment_type>,
                                    seqan3::align_cfg::mode<seqan3::detail::local_alignment_type>,
                                    seqan3::align_cfg::parallel,
                                    seqan3::align_cfg::result<>,
                                    seqan3::align_cfg::scoring<seqan3::nucleotide_scoring_scheme<int8_t>>,
                                    seqan3::detail::vectorise_tag,
                                    seqan3::detail::alignment_result_capture_element<alignment_result_t>>;

TYPED_TEST_SUITE(alignment_configuration_test, test_types, );

// TODO: this should go to a typed configuration test that also checks the alignment configuration
TEST(alignment_configuration_test, symmetric_configuration)
{
    for (uint8_t i = 0; i < static_cast<uint8_t>(seqan3::detail::align_config_id::SIZE); ++i)
    {
        // no element can occur twice in a configuration
        EXPECT_FALSE(seqan3::detail::compatibility_table<seqan3::detail::align_config_id>[i][i])
            << "There is a TRUE value on the diagonal of the search configuration matrix.";
        for (uint8_t j = 0; j < i; ++j)
        {
            // symmetric matrix
            EXPECT_EQ(seqan3::detail::compatibility_table<seqan3::detail::align_config_id>[i][j],
                      seqan3::detail::compatibility_table<seqan3::detail::align_config_id>[j][i])
                << "Search configuration matrix is not symmetric.";
        }
    }
}

TEST(alignment_configuration_test, number_of_configs)
{
    // NOTE(rrahn): You must update this test if you add a new value to seqan3::align_cfg::id
    EXPECT_EQ(static_cast<uint8_t>(seqan3::detail::align_config_id::SIZE), 12);
}

TYPED_TEST(alignment_configuration_test, config_element)
{
    EXPECT_TRUE((seqan3::detail::config_element<TypeParam>));
}

TYPED_TEST(alignment_configuration_test, configuration_exists)
{
    seqan3::configuration cfg{TypeParam{}};
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
    seqan3::configuration cfg{TypeParam{}};
    helper_exists<TypeParam, decltype(cfg)>();
}
