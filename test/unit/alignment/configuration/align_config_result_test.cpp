// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/configuration/align_config_result.hpp>
#include <seqan3/core/algorithm/configuration.hpp>

using namespace seqan3;

template <typename test_t>
struct align_cfg_result_test : public ::testing::Test
{};

using test_types = ::testing::Types<detail::with_score_type,
                                    detail::with_back_coordinate_type,
                                    detail::with_front_coordinate_type,
                                    detail::with_alignment_type>;

TYPED_TEST_SUITE(align_cfg_result_test, test_types, );

TEST(align_config_max_error, config_element)
{
    EXPECT_TRUE((detail::config_element<align_cfg::result<detail::with_score_type>>));
}

template <typename type>
auto type_to_variable()
{
    using namespace seqan3::align_cfg;

    if constexpr (std::is_same_v<type, detail::with_score_type>)
    {
        return with_score;
    }
    else if constexpr (std::is_same_v<type, detail::with_back_coordinate_type>)
    {
        return with_back_coordinate;
    }
    else if constexpr (std::is_same_v<type, detail::with_front_coordinate_type>)
    {
        return with_front_coordinate;
    }
    else
    {
        return with_alignment;
    }
}

TYPED_TEST(align_cfg_result_test, configuration)
{
    {
        align_cfg::result elem{TypeParam{}};
        configuration cfg{elem};
        EXPECT_TRUE((std::is_same_v<std::remove_reference_t<decltype(get<align_cfg::result>(cfg).value)>,
                                    TypeParam>));
    }

    {
        configuration cfg{align_cfg::result{type_to_variable<TypeParam>()}};
        EXPECT_TRUE((std::is_same_v<std::remove_reference_t<decltype(get<align_cfg::result>(cfg).value)>,
                                    TypeParam>));
    }
}

TYPED_TEST(align_cfg_result_test, score_type)
{
    // second template argument defaulted
    EXPECT_TRUE((std::is_same_v<decltype(align_cfg::result{TypeParam{}}),
                                align_cfg::result<TypeParam, int32_t>>));

    // second template argument deduced in construction
    EXPECT_TRUE((std::is_same_v<decltype(align_cfg::result{TypeParam{}, using_score_type<double>}),
                                align_cfg::result<TypeParam, double>>));

    // member type variable `score_type`
    EXPECT_TRUE((std::is_same_v<typename decltype(align_cfg::result{TypeParam{}, using_score_type<double>})::score_type,
                                double>));
}
