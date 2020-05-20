// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/detail/pack_algorithm.hpp>
#include <seqan3/search/configuration/hit.hpp>

#include "../../core/algorithm/pipeable_config_element_test_template.hpp"
#include "../../core/algorithm/pipeable_config_element_selector_test_template.hpp"

// ---------------------------------------------------------------------------------------------------------------------
// test template : pipeable_config_element_test, config_selector_test
// ---------------------------------------------------------------------------------------------------------------------

using test_types = ::testing::Types<seqan3::detail::hit_all_tag,
                                    seqan3::detail::hit_all_best_tag,
                                    seqan3::detail::hit_single_best_tag,
                                    seqan3::search_cfg::hit_strata>;

INSTANTIATE_TYPED_TEST_SUITE_P(mode_elements, pipeable_config_element_test, test_types, );

template <>
struct config_selector_test<seqan3::search_cfg::hit> : public ::testing::Test
{
    using selectee_list = seqan3::type_list<seqan3::detail::hit_all_tag,
                                            seqan3::detail::hit_all_best_tag,
                                            seqan3::detail::hit_single_best_tag,
                                            seqan3::search_cfg::hit_strata>;
};

INSTANTIATE_TYPED_TEST_SUITE_P(hit, config_selector_test, seqan3::search_cfg::hit, );

// ---------------------------------------------------------------------------------------------------------------------
// individual tests
// ---------------------------------------------------------------------------------------------------------------------

TEST(config_element_test, tags)
{
    seqan3::configuration elem_all = seqan3::search_cfg::hit_all;
    EXPECT_TRUE((std::same_as<decltype(elem_all), seqan3::configuration<seqan3::detail::hit_all_tag>>));

    seqan3::configuration elem_all_best = seqan3::search_cfg::hit_all_best;
    EXPECT_TRUE((std::same_as<decltype(elem_all_best), seqan3::configuration<seqan3::detail::hit_all_best_tag>>));

    seqan3::configuration elem_single_best = seqan3::search_cfg::hit_single_best;
    EXPECT_TRUE((std::same_as<decltype(elem_single_best), seqan3::configuration<seqan3::detail::hit_single_best_tag>>));
}

TEST(hit_strata_test, member_variable)
{
    {   // default construction
        seqan3::search_cfg::hit_strata strata_mode{};
        EXPECT_EQ(strata_mode.value, 0);
    }

    {   // construct with value
        seqan3::search_cfg::hit_strata strata_mode{3};
        EXPECT_EQ(strata_mode.value, 3);
    }

    {   // assign value
        seqan3::search_cfg::hit_strata strata_mode{};
        strata_mode.value = 3;
        EXPECT_EQ(strata_mode.value, 3);
    }
}
