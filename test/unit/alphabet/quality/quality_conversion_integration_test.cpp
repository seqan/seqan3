// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <numeric>
#include <sstream>
#include <tuple>
#include <utility>

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/all.hpp>

using namespace seqan3;

template <typename T>
using quality_conversion = ::testing::Test;

// add all alphabets from the quality sub module here
using quality_conversion_types = type_list<phred42, phred63, phred68legacy>;
using quality_conversion_gtest_types = detail::transfer_template_args_onto_t<quality_conversion_types, ::testing::Types>;

TYPED_TEST_CASE(quality_conversion, quality_conversion_gtest_types);

TYPED_TEST(quality_conversion, explicit_conversion)
{
    meta::for_each(quality_conversion_types{}, [&] (auto && qual) constexpr
    {
        using out_type = std::decay_t<decltype(qual)>;
        EXPECT_EQ(static_cast<out_type>(TypeParam{ 0}), out_type{ 0});
        EXPECT_EQ(static_cast<out_type>(TypeParam{ 5}), out_type{ 5});
        EXPECT_EQ(static_cast<out_type>(TypeParam{15}), out_type{15});
        EXPECT_EQ(static_cast<out_type>(TypeParam{20}), out_type{20});
        EXPECT_EQ(static_cast<out_type>(TypeParam{40}), out_type{40});
    });
}
