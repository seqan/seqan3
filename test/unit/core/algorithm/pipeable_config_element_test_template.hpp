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

// The function is given a list of pairs with the configuration element and the associated tabu list.
// This list is only deefined once per test case. This function allows to extract a list containing only the
// configuration element types without their associated tabu lists. This list is then used to
// instantiate the typed tests.
template <typename target_list_t, typename config_tabu_pair_list_t>
constexpr auto generate_config_element_list()
{
    if constexpr (seqan3::list_traits::size<config_tabu_pair_list_t> == 0u)
    {
        return target_list_t{};
    }
    else
    {
        using config_element_t = typename seqan3::list_traits::front<config_tabu_pair_list_t>::first_type;
        return generate_config_element_list<seqan3::list_traits::concat<target_list_t,
                                                                         seqan3::type_list<config_element_t>>,
                                            seqan3::list_traits::drop_front<config_tabu_pair_list_t>>();
    }
}

// Defines a type list with the tabu list information stripped away.
template <typename config_tabu_pair_list_t>
using pure_config_type_list = decltype(generate_config_element_list<seqan3::type_list<>,
                                                                    config_tabu_pair_list_t>());

// Generates a seqan3::configuration object that does not contain the configuration elements given
// by the tabu_list.
// The function recursively steps through the tabu list and removes the element from the config_list_t.
// The first element of the tabu list is cut of before every recursive call.
// If the tabu list is empty, a configuration object with only valid configuration elements will be returned.
template <typename config_list_t, typename tabu_list_t>
constexpr auto make_configuration()
{
    if constexpr (seqan3::list_traits::size<tabu_list_t> == 0u)
    {
        return seqan3::detail::transfer_template_args_onto_t<config_list_t, seqan3::configuration>{};
    }
    else
    {
        using tabu_config_t = seqan3::list_traits::front<tabu_list_t>;
        constexpr std::ptrdiff_t pos = seqan3::list_traits::find<tabu_config_t, config_list_t>;
        using split_list_t = seqan3::list_traits::split_after<pos, config_list_t>;
        using new_config_list_t =
            seqan3::list_traits::concat<typename split_list_t::first_type,
                                        seqan3::list_traits::drop_front<typename split_list_t::second_type>>;

        return make_configuration<new_config_list_t, seqan3::list_traits::drop_front<tabu_list_t>>();
    }
}

template <typename config_list_t, typename tabu_list_t>
using make_pipeable_configuration = decltype(make_configuration<config_list_t, tabu_list_t>());

template <typename config_fixture_t>
struct pipeable_config_element_test : public ::testing::Test, public config_fixture_t
{};

TYPED_TEST_SUITE_P(pipeable_config_element_test);

TYPED_TEST_P(pipeable_config_element_test, concept_check)
{
    EXPECT_TRUE((seqan3::config_element_specialisation<typename TestFixture::config_type>));
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
    EXPECT_TRUE((std::same_as<decltype(cfg), seqan3::configuration<config_t>>));
}

TYPED_TEST_P(pipeable_config_element_test, configuration_assignment)
{
    using config_t = typename TestFixture::config_type;

    seqan3::configuration cfg = config_t{};
    EXPECT_TRUE((std::same_as<decltype(cfg), seqan3::configuration<config_t>>));
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

template <template <typename ...> typename t, typename cfg_t>
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

    // Test combineability with every config that can be combined with.
    seqan3::detail::for_each<configs_list_t>([&](auto other_config)
    {
        using other_config_t = typename decltype(other_config)::type;

        EXPECT_TRUE((seqan3::config_element_pipeable_with<typename TestFixture::config_type, other_config_t>));
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
        using expected_config_types_t =
            seqan3::list_traits::concat<seqan3::detail::transfer_template_args_onto_t<configuration_t,
                                                                                      seqan3::type_list>,
                                        seqan3::type_list<config_t>>;
        using expected_configuration_t =
            seqan3::detail::transfer_template_args_onto_t<expected_config_types_t, seqan3::configuration>;

        {   // lvalue | lvalue
            auto cfg = compatible_configuration | elem;
            EXPECT_TRUE((std::is_same_v<decltype(cfg), expected_configuration_t>));
        }

        {   // lvalue | rvalue
            auto cfg = compatible_configuration | config_t{};
            EXPECT_TRUE((std::is_same_v<decltype(cfg), expected_configuration_t>));
        }

        {   // rvalue | lvalue
            auto cfg = configuration_t{} | elem;
            EXPECT_TRUE((std::is_same_v<decltype(cfg), expected_configuration_t>));
        }

        {   // rvalue | rvalue
            auto cfg = configuration_t{} | config_t{};
            EXPECT_TRUE((std::is_same_v<decltype(cfg), expected_configuration_t>));
        }
    }

    { // Test with config element on the left hand side.
            using expected_config_types_t =
                    seqan3::list_traits::concat<seqan3::type_list<config_t>,
                                                seqan3::detail::transfer_template_args_onto_t<configuration_t,
                                                                                              seqan3::type_list>>;
            using expected_configuration_t =
                    seqan3::detail::transfer_template_args_onto_t<expected_config_types_t, seqan3::configuration>;

        {   // lvalue | lvalue
            auto cfg = elem | compatible_configuration;
            EXPECT_TRUE((std::is_same_v<decltype(cfg), expected_configuration_t>));
        }

        {   // rvalue | lvalue
            auto cfg = config_t{} | compatible_configuration;
            EXPECT_TRUE((std::is_same_v<decltype(cfg), expected_configuration_t>));
        }

        {   // lvalue | rvalue
            auto cfg = elem | configuration_t{};
            EXPECT_TRUE((std::is_same_v<decltype(cfg), expected_configuration_t>));
        }

        {   // rvalue | rvalue
            auto cfg = config_t{} | configuration_t{};
            EXPECT_TRUE((std::is_same_v<decltype(cfg), expected_configuration_t>));
        }
    }

    { // Test with empty configuration.
        using expected_configuration_t = seqan3::configuration<config_t>;

        {   // lvalue | empty
            auto cfg = elem | seqan3::configuration<>{};
            EXPECT_TRUE((std::is_same_v<decltype(cfg), expected_configuration_t>));
        }

        {   // rvalue | empty
            auto cfg = config_t{} | seqan3::configuration<>{};
            EXPECT_TRUE((std::is_same_v<decltype(cfg), expected_configuration_t>));
        }

        {   // empty | lvalue
            auto cfg = seqan3::configuration<>{} | elem;
            EXPECT_TRUE((std::is_same_v<decltype(cfg), expected_configuration_t>));
        }

        {   // empty | rvalue
            auto cfg = seqan3::configuration<>{} | config_t{};
            EXPECT_TRUE((std::is_same_v<decltype(cfg), expected_configuration_t>));
        }
    }
}

TYPED_TEST_P(pipeable_config_element_test, invalid_pipeability)
{
    using config_t = typename TestFixture::config_type;
    using incompatible_config_same_domain_t = seqan3::list_traits::front<typename TestFixture::tabu_list_type>;

    EXPECT_FALSE((seqan3::config_element_pipeable_with<config_t, incompatible_config_same_domain_t>));
    EXPECT_FALSE((seqan3::config_element_pipeable_with<incompatible_config_same_domain_t, config_t>));
    EXPECT_FALSE((seqan3::config_element_pipeable_with<config_t, foo>));
    EXPECT_FALSE((seqan3::config_element_pipeable_with<foo, config_t>));
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
