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

using namespace seqan3;

template <typename type>
struct align_cfg_mode_test : public ::testing::Test
{};

using test_types = ::testing::Types<detail::global_alignment_type, detail::local_alignment_type>;
TYPED_TEST_SUITE(align_cfg_mode_test, test_types, );

TYPED_TEST(align_cfg_mode_test, config_element)
{
    EXPECT_EQ(detail::config_element<align_cfg::mode<TypeParam>>, true);
}

TYPED_TEST(align_cfg_mode_test, configuration)
{
    {
        align_cfg::mode elem{TypeParam{}};
        configuration cfg(elem);
        EXPECT_EQ((std::is_same_v<std::remove_reference_t<decltype(get<align_cfg::mode>(cfg).value)>,
                                  TypeParam>), true);
    }

    {
        configuration cfg{align_cfg::mode{TypeParam{}}};
        EXPECT_EQ((std::is_same_v<std::remove_reference_t<decltype(get<align_cfg::mode>(cfg).value)>,
                                  TypeParam>), true);
    }
}

template <typename type>
constexpr auto get_inline_variable()
{
    if constexpr (std::is_same_v<type, detail::global_alignment_type>)
    {
        return global_alignment;
    }
    else if constexpr (std::is_same_v<type, detail::local_alignment_type>)
    {
        return local_alignment;
    }
    else
    {
        return std::ignore;
    }
}

TYPED_TEST(align_cfg_mode_test, construction_from_variable)
{
    configuration cfg{align_cfg::mode{get_inline_variable<TypeParam>()}};
    EXPECT_EQ((std::is_same_v<std::remove_reference_t<decltype(get<align_cfg::mode>(cfg).value)>,
                              TypeParam>), true);
}
