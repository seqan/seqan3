// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

// This file can be used to test the basic functionality of a configuration element selector

#include <gtest/gtest.h>

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/detail/pack_algorithm.hpp>

template <typename selector_t>
struct config_selector_test : public ::testing::Test
{};

TYPED_TEST_SUITE_P(config_selector_test);

TYPED_TEST_P(config_selector_test, concept_check)
{
    EXPECT_TRUE((seqan3::detail::config_element<TypeParam>));
    foo<TypeParam>();
}

TYPED_TEST_P(config_selector_test, construction_from_element)
{
    seqan3::detail::for_each<typename TestFixture::selectee_list>([&] (auto selected_identity_type)
    {
        using selected_type = typename decltype(selected_identity_type)::type; // remove type_identity wrapper

        TypeParam{selected_type{}};
    });
}

TYPED_TEST_P(config_selector_test, asignment_from_element)
{
    seqan3::detail::for_each<typename TestFixture::selectee_list>([&] (auto selected_identity_type)
    {
        using selected_type = typename decltype(selected_identity_type)::type; // remove type_identity wrapper

        TypeParam select{};
        select = selected_type{};
    });
}

namespace seqan3::detail
{
struct test_accessor
{
    template <typename expected_type, typename selector_type>
    static void check_selection(selector_type const & selector)
    {
        std::visit([] (auto selected)
        {
            EXPECT_TRUE((std::same_as<expected_type, decltype(selected)>));
        }, selector.selection);
    }
};
} // seqan3::detail

TYPED_TEST_P(config_selector_test, member_variable_selection)
{
    seqan3::detail::for_each<typename TestFixture::selectee_list>([&] (auto selected_identity_type)
    {
        using selected_type = typename decltype(selected_identity_type)::type; // remove type_identity wrapper

        TypeParam m{selected_type{}};

        seqan3::detail::test_accessor::check_selection<selected_type>(m);
    });
}

REGISTER_TYPED_TEST_SUITE_P(config_selector_test,
                            concept_check,
                            construction_from_element,
                            asignment_from_element,
                            member_variable_selection);
