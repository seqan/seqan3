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

template <typename test_t>
struct align_cfg_result_test : public ::testing::Test
{};

using test_types = ::testing::Types<seqan3::detail::with_score_type,
                                    seqan3::detail::with_back_coordinate_type,
                                    seqan3::detail::with_front_coordinate_type,
                                    seqan3::detail::with_alignment_type>;

TYPED_TEST_SUITE(align_cfg_result_test, test_types, );

TEST(align_config_max_error, config_element)
{
    EXPECT_TRUE((seqan3::detail::config_element<seqan3::align_cfg::result<seqan3::detail::with_score_type>>));
}

template <typename type>
auto type_to_variable()
{
    if constexpr (std::is_same_v<type, seqan3::detail::with_score_type>)
    {
        return seqan3::with_score;
    }
    else if constexpr (std::is_same_v<type, seqan3::detail::with_back_coordinate_type>)
    {
        return seqan3::with_back_coordinate;
    }
    else if constexpr (std::is_same_v<type, seqan3::detail::with_front_coordinate_type>)
    {
        return seqan3::with_front_coordinate;
    }
    else
    {
        return seqan3::with_alignment;
    }
}

TYPED_TEST(align_cfg_result_test, configuration)
{
    {
        seqan3::align_cfg::result elem{TypeParam{}};
        seqan3::configuration cfg{elem};
        EXPECT_TRUE((std::is_same_v<std::remove_reference_t<decltype(seqan3::get<seqan3::align_cfg::result>(cfg).value)>,
                                    TypeParam>));
    }

    {
        seqan3::configuration cfg{seqan3::align_cfg::result{type_to_variable<TypeParam>()}};
        EXPECT_TRUE((std::is_same_v<std::remove_reference_t<decltype(seqan3::get<seqan3::align_cfg::result>(cfg).value)>,
                                    TypeParam>));
    }
}

TYPED_TEST(align_cfg_result_test, score_type)
{
    using type_of_score = decltype(seqan3::align_cfg::result{TypeParam{}, seqan3::using_score_type<double>});

    // second template argument defaulted
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::align_cfg::result{TypeParam{}}),
                 seqan3::align_cfg::result<TypeParam, int32_t>>));

    // second template argument deduced in construction
    EXPECT_TRUE((std::is_same_v<type_of_score, seqan3::align_cfg::result<TypeParam, double>>));

    // member type variable `score_type`
    EXPECT_TRUE((std::is_same_v<typename type_of_score::score_type, double>));
}
