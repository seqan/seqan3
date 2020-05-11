// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_mode.hpp>
#include <seqan3/core/algorithm/configuration.hpp>

#include "../../core/algorithm/pipeable_config_element_test_template.hpp"

// ---------------------------------------------------------------------------------------------------------------------
// pipeable_config_element_test template
// ---------------------------------------------------------------------------------------------------------------------

using config_element_types = ::testing::Types<seqan3::align_cfg::mode<seqan3::detail::global_alignment_type>,
                                              seqan3::align_cfg::mode<seqan3::detail::local_alignment_type>>;

INSTANTIATE_TYPED_TEST_SUITE_P(mode, pipeable_config_element_test, config_element_types, );

// ---------------------------------------------------------------------------------------------------------------------
// align_cfg_mode_test
// ---------------------------------------------------------------------------------------------------------------------

template <typename type>
struct align_cfg_mode_test : public ::testing::Test
{};

using test_types = ::testing::Types<seqan3::detail::global_alignment_type, seqan3::detail::local_alignment_type>;
TYPED_TEST_SUITE(align_cfg_mode_test, test_types, );

TYPED_TEST(align_cfg_mode_test, configuration)
{
    {
        seqan3::align_cfg::mode elem{TypeParam{}};
        seqan3::configuration cfg(elem);
        EXPECT_EQ((std::is_same_v<std::remove_reference_t<decltype(seqan3::get<seqan3::align_cfg::mode>(cfg).value)>,
                                  TypeParam>), true);
    }

    {
        seqan3::configuration cfg{seqan3::align_cfg::mode{TypeParam{}}};
        EXPECT_EQ((std::is_same_v<std::remove_reference_t<decltype(seqan3::get<seqan3::align_cfg::mode>(cfg).value)>,
                                  TypeParam>), true);
    }
}

template <typename type>
constexpr auto get_inline_variable()
{
    if constexpr (std::is_same_v<type, seqan3::detail::global_alignment_type>)
    {
        return seqan3::global_alignment;
    }
    else if constexpr (std::is_same_v<type, seqan3::detail::local_alignment_type>)
    {
        return seqan3::local_alignment;
    }
    else
    {
        return std::ignore;
    }
}

TYPED_TEST(align_cfg_mode_test, construction_from_variable)
{
    seqan3::configuration cfg{seqan3::align_cfg::mode{get_inline_variable<TypeParam>()}};
    EXPECT_EQ((std::is_same_v<std::remove_reference_t<decltype(seqan3::get<seqan3::align_cfg::mode>(cfg).value)>,
                              TypeParam>), true);
}
