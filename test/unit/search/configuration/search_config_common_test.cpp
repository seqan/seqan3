// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/search/configuration/all.hpp>
#include <seqan3/search/search_result.hpp>
#include <seqan3/utility/type_list/traits.hpp>

#include "../../core/configuration/pipeable_config_element_test_template.hpp"

// Define some aliases to make the list below more readable.
namespace cfg = seqan3::search_cfg;

using nil_t = seqan3::detail::empty_type;
using search_result_t = seqan3::search_result<nil_t, nil_t, nil_t, nil_t>;
using all_hit_configs_t =
    seqan3::type_list<cfg::hit, cfg::hit_all, cfg::hit_all_best, cfg::hit_single_best, cfg::hit_strata>;
// Needed for the on result config
[[maybe_unused]] auto on_result_callback = []([[maybe_unused]] auto && res) {};
using callback_t = decltype(on_result_callback);

// A list of config types to test, associated with their incompatible config classes defined as a taboo list.
// We later use this taboo list to generate a configuration object containing only the valid combinations for each
// test type.
using search_config_and_taboo_types = seqan3::type_list<
    // hit configs
    std::pair<cfg::hit, all_hit_configs_t>,
    std::pair<cfg::hit_all, all_hit_configs_t>,
    std::pair<cfg::hit_all_best, all_hit_configs_t>,
    std::pair<cfg::hit_single_best, all_hit_configs_t>,
    std::pair<cfg::hit_strata, all_hit_configs_t>,
    // max error configs
    std::pair<cfg::max_error_total, seqan3::type_list<cfg::max_error_total>>,
    std::pair<cfg::max_error_substitution, seqan3::type_list<cfg::max_error_substitution>>,
    std::pair<cfg::max_error_deletion, seqan3::type_list<cfg::max_error_deletion>>,
    std::pair<cfg::max_error_insertion, seqan3::type_list<cfg::max_error_insertion>>,
    // output configs
    std::pair<cfg::output_query_id, seqan3::type_list<cfg::output_query_id>>,
    std::pair<cfg::output_reference_id, seqan3::type_list<cfg::output_reference_id>>,
    std::pair<cfg::output_reference_begin_position, seqan3::type_list<cfg::output_reference_begin_position>>,
    std::pair<cfg::output_index_cursor, seqan3::type_list<cfg::output_index_cursor>>,
    // other configs
    std::pair<cfg::parallel, seqan3::type_list<cfg::parallel>>,
    std::pair<cfg::on_result<callback_t>, seqan3::type_list<cfg::on_result<callback_t>>>,
    std::pair<cfg::detail::result_type<search_result_t>, seqan3::type_list<cfg::detail::result_type<search_result_t>>>>;

// The pure list of configuration elements to instantiate the typed test case with.
using search_config_types = pure_config_type_list<search_config_and_taboo_types>;

template <typename config_t>
class test_fixture
{
public:
    // The actual config type that is tested.
    using config_type = config_t;

    // The taboo list associated with the given TypeParam element.
    using taboo_list_type = typename seqan3::list_traits::at<
        seqan3::list_traits::find<config_t, search_config_types>, // determine the index
        search_config_and_taboo_types>::second_type;              // extract the taboo list.

    // A compatible configuration type for the current configuration element to test.
    using compatible_configuration_type = make_pipeable_configuration<search_config_types, taboo_list_type>;

    // The type of the configuration element ids.
    using config_id_type = seqan3::detail::search_config_id;

    // NOTE: You must update this number if you add a new entity to seqan3::detail::search_config_id.
    // config_count is used to check that the config size is correct.
    // And don't forget to add the new config into the above test fixture (via search_config_and_taboo_types).
    static constexpr int8_t config_count = 12;
};

// Configuration element type list as gtest suitable testing::Types
using fixture_types = seqan3::list_traits::transform<test_fixture, search_config_types>;
using test_types = seqan3::detail::transfer_template_args_onto_t<fixture_types, ::testing::Types>;

INSTANTIATE_TYPED_TEST_SUITE_P(search_configuration_test, pipeable_config_element_test, test_types, );
