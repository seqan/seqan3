// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

// This file can be used to test the basic functionality of a configuration element

#include <gtest/gtest.h>

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

#include "configuration_mock.hpp"

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
    EXPECT_TRUE((std::same_as<decltype(cfg), seqan3::configuration<TypeParam>>));
}

TYPED_TEST_P(pipeable_config_element_test, configuration_assignment)
{
    seqan3::configuration cfg = TypeParam{};
    EXPECT_TRUE((std::same_as<decltype(cfg), seqan3::configuration<TypeParam>>));
}

// make mock config element foo always pipeable with any other config element
namespace seqan3::detail
{

template <typename t>
struct is_configuration_valid<t, foo> : public std::true_type
{};

template <typename t>
struct is_configuration_valid<foo, t> : public std::true_type
{};

}

TYPED_TEST_P(pipeable_config_element_test, pipeability)
{
    TypeParam elem{};
    foo dummy{};

    {   // lvalue | lvalue
        auto cfg = dummy | elem;
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<foo, TypeParam>>));
    }

    {   // rvalue | lvalue
        auto cfg = TypeParam{} | dummy;
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<TypeParam, foo>>));
    }

    {   // lvalue | rvalue
        auto cfg = dummy | TypeParam{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<foo, TypeParam>>));
    }

    {   // rvalue | rvalue
        auto cfg = foo{} | TypeParam{};
        EXPECT_TRUE((std::is_same_v<decltype(cfg), seqan3::configuration<foo, TypeParam>>));
    }
}

REGISTER_TYPED_TEST_SUITE_P(pipeable_config_element_test,
                            concept_check,
                            standard_construction,
                            configuration_construction,
                            configuration_assignment,
                            pipeability);
