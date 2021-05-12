// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

// This file can be used to test the basic functionality of a configuration element

#pragma once

#include <gtest/gtest.h>

#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>
#include <seqan3/test/expect_same_type.hpp>

#include "../configuration/configuration_mock.hpp"

template <typename format_t>
struct pipeable_config_element_test : public ::testing::Test
{};

TYPED_TEST_SUITE_P(pipeable_config_element_test);

TYPED_TEST_P(pipeable_config_element_test, concept_check)
{
    EXPECT_TRUE((seqan3::detail::config_element<TypeParam>));
}

TYPED_TEST_P(pipeable_config_element_test, standard_construction)
{
    EXPECT_TRUE((std::is_default_constructible_v<TypeParam>));
    EXPECT_TRUE((std::is_copy_constructible_v<TypeParam>));
    EXPECT_TRUE((std::is_move_constructible_v<TypeParam>));
    EXPECT_TRUE((std::is_copy_assignable_v<TypeParam>));
    EXPECT_TRUE((std::is_move_assignable_v<TypeParam>));
}

TYPED_TEST_P(pipeable_config_element_test, configuration_construction)
{
    seqan3::configuration cfg{TypeParam{}};
    EXPECT_SAME_TYPE(decltype(cfg), seqan3::configuration<TypeParam>);
}

TYPED_TEST_P(pipeable_config_element_test, configuration_assignment)
{
    seqan3::configuration cfg = TypeParam{};
    EXPECT_SAME_TYPE(decltype(cfg), seqan3::configuration<TypeParam>);
}

TYPED_TEST_P(pipeable_config_element_test, exists)
{
    using configuration_t = seqan3::configuration<TypeParam>;

    EXPECT_TRUE(configuration_t::template exists<TypeParam>());
}

REGISTER_TYPED_TEST_SUITE_P(pipeable_config_element_test,
                            concept_check,
                            standard_construction,
                            configuration_construction,
                            configuration_assignment,
                            exists);
