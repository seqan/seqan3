// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <numeric>
#include <sstream>
#include <tuple>
#include <utility>

#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/phred63.hpp>
#include <seqan3/alphabet/quality/phred68solexa.hpp>
#include <seqan3/alphabet/quality/phred94.hpp>
#include <seqan3/utility/type_list/detail/type_list_algorithm.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

template <typename T>
using quality_conversion = ::testing::Test;

// add all alphabets from the quality sub module here
using quality_conversion_types =
    seqan3::type_list<seqan3::phred42, seqan3::phred63, seqan3::phred68solexa, seqan3::phred94>;
using quality_conversion_gtest_types =
    seqan3::detail::transfer_template_args_onto_t<quality_conversion_types, ::testing::Types>;

TYPED_TEST_SUITE(quality_conversion, quality_conversion_gtest_types, );

TYPED_TEST(quality_conversion, explicit_conversion)
{
    seqan3::detail::for_each<quality_conversion_types>(
        [&](auto qual) constexpr
        {
            using out_type = std::decay_t<typename decltype(qual)::type>;
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_phred(0)), out_type{}.assign_phred(0));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_phred(5)), out_type{}.assign_phred(5));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_phred(15)), out_type{}.assign_phred(15));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_phred(20)), out_type{}.assign_phred(20));
            EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_phred(40)), out_type{}.assign_phred(40));
        });
}
