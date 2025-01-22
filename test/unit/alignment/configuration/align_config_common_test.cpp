// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_band.hpp>
#include <seqan3/alignment/configuration/align_config_debug.hpp>
#include <seqan3/alignment/configuration/align_config_gap_cost_affine.hpp>
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_min_score.hpp>
#include <seqan3/alignment/configuration/align_config_on_result.hpp>
#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/alignment/configuration/align_config_parallel.hpp>
#include <seqan3/alignment/configuration/align_config_result_type.hpp>
#include <seqan3/alignment/configuration/align_config_score_type.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/configuration/align_config_vectorised.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/utility/type_list/traits.hpp>

#include "../../core/configuration/pipeable_config_element_test_template.hpp"

// Define some aliases to make the list below more readable.
namespace cfg = seqan3::align_cfg;

using alignment_result_t = seqan3::alignment_result<seqan3::detail::alignment_result_value_type<int, int, int>>;
using nt_scheme = seqan3::nucleotide_scoring_scheme<int8_t>;

// Needed for the on result config
[[maybe_unused]] auto on_result_callback = []([[maybe_unused]] auto && res) {};
using callback_t = decltype(on_result_callback);

// A list of config types to test, associated with their incompatible config classes defined as a taboo list.
// We later use this taboo list to generate a configuration object containing only the valid combinations for each
// test type.
using align_config_and_taboo_types = seqan3::type_list<
    // method configs
    std::pair<cfg::method_global, seqan3::type_list<cfg::method_global, cfg::method_local>>,
    std::pair<cfg::method_local, seqan3::type_list<cfg::method_local, cfg::method_global, cfg::min_score>>,
    // output configs
    std::pair<cfg::output_sequence1_id, seqan3::type_list<cfg::output_sequence1_id>>,
    std::pair<cfg::output_sequence2_id, seqan3::type_list<cfg::output_sequence2_id>>,
    std::pair<cfg::output_score, seqan3::type_list<cfg::output_score>>,
    std::pair<cfg::output_begin_position, seqan3::type_list<cfg::output_begin_position>>,
    std::pair<cfg::output_end_position, seqan3::type_list<cfg::output_end_position>>,
    std::pair<cfg::output_alignment, seqan3::type_list<cfg::output_alignment>>,
    // other configs
    std::pair<cfg::band_fixed_size, seqan3::type_list<cfg::band_fixed_size>>,
    std::pair<cfg::detail::debug, seqan3::type_list<cfg::detail::debug>>,
    std::pair<cfg::gap_cost_affine, seqan3::type_list<cfg::gap_cost_affine>>,
    std::pair<cfg::min_score, seqan3::type_list<cfg::min_score, cfg::method_local>>,
    std::pair<cfg::on_result<callback_t>, seqan3::type_list<cfg::on_result<callback_t>>>,
    std::pair<cfg::parallel, seqan3::type_list<cfg::parallel>>,
    std::pair<cfg::detail::result_type<alignment_result_t>,
              seqan3::type_list<cfg::detail::result_type<alignment_result_t>>>,
    std::pair<cfg::score_type<int32_t>, seqan3::type_list<cfg::score_type<int32_t>>>,
    std::pair<cfg::scoring_scheme<nt_scheme>, seqan3::type_list<cfg::scoring_scheme<nt_scheme>>>,
    std::pair<cfg::vectorised, seqan3::type_list<cfg::vectorised>>>;

// The pure list of configuration elements to instantiate the typed test case with.
using align_config_types = pure_config_type_list<align_config_and_taboo_types>;

template <typename config_t>
class test_fixture
{
public:
    // The actual config type that is tested.
    using config_type = config_t;
    // The taboo list associated with the given TypeParam element.
    using taboo_list_type = typename seqan3::list_traits::at<
        seqan3::list_traits::find<config_type, align_config_types>, // determine the index
        align_config_and_taboo_types>::second_type;                 // extract the taboo list.

    // A compatible configuration type for the current configuration element to test.
    using compatible_configuration_type = make_pipeable_configuration<align_config_types, taboo_list_type>;

    using config_id_type = seqan3::detail::align_config_id;

    // NOTE: You must update this number if you add a new entity to seqan3::detail::align_config_id.
    // config_count is used to check that the config size is correct.
    // And don't forget to add the new config into the above test fixture (via align_config_and_taboo_types).
    static constexpr int8_t config_count = 18;
};

// Configuration element type list as gtest suitable testing::Types
using fixture_types = seqan3::list_traits::transform<test_fixture, align_config_types>;
using test_types = seqan3::detail::transfer_template_args_onto_t<fixture_types, ::testing::Types>;

TYPED_TEST_SUITE(alignment_configuration_fixture, test_types, );

INSTANTIATE_TYPED_TEST_SUITE_P(alignment_configuration_test, pipeable_config_element_test, test_types, );
