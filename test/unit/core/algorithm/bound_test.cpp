// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/core/algorithm/bound.hpp>

template <typename t>
struct bound_test : ::testing::Test
{};

using testing_types = ::testing::Types<int8_t, int16_t, uint32_t, float>;
TYPED_TEST_SUITE(bound_test, testing_types, );

TYPED_TEST(bound_test, lower_bound)
{
    seqan3::lower_bound lb{static_cast<TypeParam>(-5)};
    EXPECT_TRUE((std::is_same_v<decltype(lb), seqan3::lower_bound<TypeParam>>));
    EXPECT_NEAR(lb.get(), static_cast<TypeParam>(-5), 0.1);
}

TYPED_TEST(bound_test, upper_bound)
{
    seqan3::upper_bound lb{static_cast<TypeParam>(5)};
    EXPECT_TRUE((std::is_same_v<decltype(lb), seqan3::upper_bound<TypeParam>>));
    EXPECT_NEAR(lb.get(), static_cast<TypeParam>(5), 0.1);
}
