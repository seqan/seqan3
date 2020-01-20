// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/search/fm_index/all.hpp>

#include <gtest/gtest.h>

using namespace sdsl;

template <typename T>
class sdsl_index_test : public ::testing::Test
{};

template <typename alphabet_strategy_t>
using sdsl_index = csa_wt<wt_blcd<bit_vector, rank_support_v<>, select_support_scan<>, select_support_scan<0>>,
                          16,
                          10000000,
                          sa_order_sa_sampling<>,
                          isa_sampling<>,
                          alphabet_strategy_t>;

using alphabet_strategy_types = ::testing::Types<byte_alphabet, plain_byte_alphabet>;

TYPED_TEST_SUITE(sdsl_index_test, alphabet_strategy_types, );

TYPED_TEST(sdsl_index_test, concepts)
{
    EXPECT_TRUE(seqan3::detail::sdsl_index<sdsl_index<TypeParam>>);
}
