// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <sstream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/detail/pack_algorithm.hpp>
#include <seqan3/core/type_list/type_list.hpp>

using namespace seqan3;

//-----------------------------------------------------------------------------
// all_of
//-----------------------------------------------------------------------------

struct is_integral_fn
{
    bool operator()(...) { return false; }

    template <typename identity_t>
        requires std::integral<typename identity_t::type>
    bool operator()(identity_t) { return true; }
};

auto is_value_type_integral = [] (auto value)
{
    return std::is_integral_v<decltype(value)>;
};

TEST(pack_algorithm, all_of_in_type_list)
{
    EXPECT_TRUE(detail::all_of<type_list<>>(is_integral_fn{}));
    EXPECT_TRUE((detail::all_of<type_list<int8_t, int16_t, uint32_t>>(is_integral_fn{})));
    EXPECT_FALSE((detail::all_of<type_list<int8_t, int16_t, uint32_t, float>>(is_integral_fn{})));
}

TEST(pack_algorithm, all_of_values)
{
    EXPECT_TRUE(detail::all_of(is_value_type_integral));
    EXPECT_TRUE((detail::all_of(is_value_type_integral, int8_t{}, int16_t{}, uint32_t{})));
    EXPECT_FALSE((detail::all_of(is_value_type_integral, int8_t{}, int16_t{}, uint32_t{}, float{})));
}
