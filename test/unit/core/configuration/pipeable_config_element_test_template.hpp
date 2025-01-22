// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

// This file can be used to test the basic functionality of a configuration element

#pragma once

#include <gtest/gtest.h>

#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/type_list/detail/type_list_algorithm.hpp>

#include "../configuration/configuration_mock.hpp"

// The function is given a list of pairs with the configuration element and the associated taboo list.
// This list is only defined once per test case. This function allows to extract a list containing only the
// configuration element types without their associated taboo lists. This list is then used to
// instantiate the typed tests.
template <typename... taboo_pairs_t>
constexpr auto pure_config_type_list_fn(seqan3::type_list<taboo_pairs_t...>)
    -> seqan3::type_list<typename taboo_pairs_t::first_type...>;

// Defines a type list with the taboo list information stripped away.
template <typename config_taboo_pair_list_t>
using pure_config_type_list = decltype(pure_config_type_list_fn(config_taboo_pair_list_t{}));

// Generates a seqan3::configuration object that does not contain the configuration elements given
// by the taboo list.
// The function recursively steps through the taboo list and removes the element from the config_list_t.
// The first element of the taboo list is cut of before every recursive call.
// If the taboo list is empty, a configuration object with only valid configuration elements will be returned.
template <typename config_list_t, typename taboo_list_t>
using make_pipeable_configuration = seqan3::detail::transfer_template_args_onto_t<
    decltype(seqan3::list_traits::detail::type_list_difference(config_list_t{}, taboo_list_t{})),
    seqan3::configuration>;

template <typename config_fixture_t>
struct pipeable_config_element_test : public ::testing::Test, public config_fixture_t
{};

TYPED_TEST_SUITE_P(pipeable_config_element_test);

TYPED_TEST_P(pipeable_config_element_test, concept_check)
{
    EXPECT_TRUE((seqan3::detail::config_element<typename TestFixture::config_type>));
}

TYPED_TEST_P(pipeable_config_element_test, standard_construction)
{
    using config_t = typename TestFixture::config_type;

    EXPECT_TRUE((std::is_default_constructible_v<config_t>));
    EXPECT_TRUE((std::is_copy_constructible_v<config_t>));
    EXPECT_TRUE((std::is_move_constructible_v<config_t>));
    EXPECT_TRUE((std::is_copy_assignable_v<config_t>));
    EXPECT_TRUE((std::is_move_assignable_v<config_t>));
}

TYPED_TEST_P(pipeable_config_element_test, configuration_construction)
{
    using config_t = typename TestFixture::config_type;

    seqan3::configuration cfg{config_t{}};
    EXPECT_SAME_TYPE(decltype(cfg), seqan3::configuration<config_t>);
}

TYPED_TEST_P(pipeable_config_element_test, configuration_assignment)
{
    using config_t = typename TestFixture::config_type;

    seqan3::configuration cfg = config_t{};
    EXPECT_SAME_TYPE(decltype(cfg), seqan3::configuration<config_t>);
}

TYPED_TEST_P(pipeable_config_element_test, symmetric_configuration)
{
    using config_id_t = typename TestFixture::config_id_type;
    for (uint8_t i = 0; i < static_cast<uint8_t>(config_id_t::SIZE); ++i)
    {
        // no element can occur twice in a configuration
        EXPECT_FALSE(seqan3::detail::compatibility_table<config_id_t>[i][i])
            << "There is a TRUE value on the diagonal of the search configuration matrix.";
        for (uint8_t j = 0; j < i; ++j)
        {
            // symmetric matrix
            EXPECT_EQ(seqan3::detail::compatibility_table<config_id_t>[i][j],
                      seqan3::detail::compatibility_table<config_id_t>[j][i])
                << "Search configuration matrix is not symmetric.";
        }
    }
}

TYPED_TEST_P(pipeable_config_element_test, number_of_configs)
{
    using config_id_t = typename TestFixture::config_id_type;

    EXPECT_EQ(static_cast<uint8_t>(config_id_t::SIZE), TestFixture::config_count);
}

TYPED_TEST_P(pipeable_config_element_test, exists)
{
    using config_t = typename TestFixture::config_type;
    using configuration_t = seqan3::configuration<config_t>;

    EXPECT_TRUE(configuration_t::template exists<config_t>());
}

template <typename t, typename cfg_t>
void helper_exists()
{
    EXPECT_TRUE(cfg_t::template exists<t>());
}

template <template <typename...> typename t, typename cfg_t>
void helper_exists()
{
    EXPECT_TRUE(cfg_t::template exists<t>());
}

TYPED_TEST_P(pipeable_config_element_test, exists_template)
{
    using config_t = typename TestFixture::config_type;

    seqan3::configuration cfg{config_t{}};
    helper_exists<config_t, decltype(cfg)>();
}

TYPED_TEST_P(pipeable_config_element_test, combineable_with)
{
    using configs_list_t =
        seqan3::detail::transfer_template_args_onto_t<typename TestFixture::compatible_configuration_type,
                                                      seqan3::type_list>;

    // Test combinability with every config that can be combined with.
    seqan3::detail::for_each<configs_list_t>(
        [&](auto other_config)
        {
            using other_config_t = typename decltype(other_config)::type;

            EXPECT_TRUE(
                (seqan3::detail::config_element_pipeable_with<typename TestFixture::config_type, other_config_t>));
        });
}

TYPED_TEST_P(pipeable_config_element_test, pipeability)
{
    using config_t = typename TestFixture::config_type;
    using configuration_t = typename TestFixture::compatible_configuration_type;

    // Test actual pipe operation.
    configuration_t compatible_configuration{};
    config_t elem{};

    { // Test with config element on the right hand side.
        using expected_config_types_t = seqan3::list_traits::concat<
            seqan3::detail::transfer_template_args_onto_t<configuration_t, seqan3::type_list>,
            seqan3::type_list<config_t>>;
        using expected_configuration_t =
            seqan3::detail::transfer_template_args_onto_t<expected_config_types_t, seqan3::configuration>;

        { // lvalue | lvalue
            auto cfg = compatible_configuration | elem;
            EXPECT_SAME_TYPE(decltype(cfg), expected_configuration_t);
        }

        { // lvalue | rvalue
            auto cfg = compatible_configuration | config_t{};
            EXPECT_SAME_TYPE(decltype(cfg), expected_configuration_t);
        }

        { // rvalue | lvalue
            auto cfg = configuration_t{} | elem;
            EXPECT_SAME_TYPE(decltype(cfg), expected_configuration_t);
        }

        { // rvalue | rvalue
            auto cfg = configuration_t{} | config_t{};
            EXPECT_SAME_TYPE(decltype(cfg), expected_configuration_t);
        }
    }

    { // Test with config element on the left hand side.
        using expected_config_types_t = seqan3::list_traits::concat<
            seqan3::type_list<config_t>,
            seqan3::detail::transfer_template_args_onto_t<configuration_t, seqan3::type_list>>;
        using expected_configuration_t =
            seqan3::detail::transfer_template_args_onto_t<expected_config_types_t, seqan3::configuration>;

        { // lvalue | lvalue
            auto cfg = elem | compatible_configuration;
            EXPECT_SAME_TYPE(decltype(cfg), expected_configuration_t);
        }

        { // rvalue | lvalue
            auto cfg = config_t{} | compatible_configuration;
            EXPECT_SAME_TYPE(decltype(cfg), expected_configuration_t);
        }

        { // lvalue | rvalue
            auto cfg = elem | configuration_t{};
            EXPECT_SAME_TYPE(decltype(cfg), expected_configuration_t);
        }

        { // rvalue | rvalue
            auto cfg = config_t{} | configuration_t{};
            EXPECT_SAME_TYPE(decltype(cfg), expected_configuration_t);
        }
    }

    { // Test with empty configuration.
        using expected_configuration_t = seqan3::configuration<config_t>;

        { // lvalue | empty
            auto cfg = elem | seqan3::configuration<>{};
            EXPECT_SAME_TYPE(decltype(cfg), expected_configuration_t);
        }

        { // rvalue | empty
            auto cfg = config_t{} | seqan3::configuration<>{};
            EXPECT_SAME_TYPE(decltype(cfg), expected_configuration_t);
        }

        { // empty | lvalue
            auto cfg = seqan3::configuration<>{} | elem;
            EXPECT_SAME_TYPE(decltype(cfg), expected_configuration_t);
        }

        { // empty | rvalue
            auto cfg = seqan3::configuration<>{} | config_t{};
            EXPECT_SAME_TYPE(decltype(cfg), expected_configuration_t);
        }
    }
}

TYPED_TEST_P(pipeable_config_element_test, invalid_pipeability)
{
    using config_t = typename TestFixture::config_type;
    using incompatible_config_same_domain_t = seqan3::list_traits::front<typename TestFixture::taboo_list_type>;

    EXPECT_FALSE((seqan3::detail::config_element_pipeable_with<config_t, incompatible_config_same_domain_t>));
    EXPECT_FALSE((seqan3::detail::config_element_pipeable_with<incompatible_config_same_domain_t, config_t>));
    EXPECT_FALSE((seqan3::detail::config_element_pipeable_with<config_t, foo>));
    EXPECT_FALSE((seqan3::detail::config_element_pipeable_with<foo, config_t>));
}

REGISTER_TYPED_TEST_SUITE_P(pipeable_config_element_test,
                            concept_check,
                            standard_construction,
                            configuration_construction,
                            configuration_assignment,
                            symmetric_configuration,
                            number_of_configs,
                            exists,
                            exists_template,
                            combineable_with,
                            pipeability,
                            invalid_pipeability);
