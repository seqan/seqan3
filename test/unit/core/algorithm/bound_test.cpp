// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/core/algorithm/bound.hpp>

using namespace seqan3;

template <typename t>
struct bound_test : ::testing::Test
{};

using testing_types = ::testing::Types<int8_t, int16_t, uint32_t, float>;
TYPED_TEST_CASE(bound_test, testing_types);

TYPED_TEST(bound_test, lower_bound)
{
    lower_bound lb{static_cast<TypeParam>(-5)};
    EXPECT_TRUE((std::is_same_v<decltype(lb), lower_bound<TypeParam>>));
    EXPECT_NEAR(lb.get(), static_cast<TypeParam>(-5), 0.1);
}

TYPED_TEST(bound_test, upper_bound)
{
    upper_bound lb{static_cast<TypeParam>(5)};
    EXPECT_TRUE((std::is_same_v<decltype(lb), upper_bound<TypeParam>>));
    EXPECT_NEAR(lb.get(), static_cast<TypeParam>(5), 0.1);
}
