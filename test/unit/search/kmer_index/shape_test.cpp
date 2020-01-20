// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/kmer_index/shape.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

using namespace seqan3;

constexpr bool construction_test()
{
    bool res{true};

    shape s1{bin_literal{0b1011}};
    res &= std::ranges::size(s1) == 4;
    res &= s1.all() == false;

    shape s2{0b1011_shape};
    res &= std::ranges::size(s2) == 4;
    res &= s2.all() == false;

    shape s3{ungapped{3}};
    res &= std::ranges::size(s3) == 3;
    res &= s3.all() == true;

    shape s4{bin_literal{0b1111}};
    res &= std::ranges::size(s4) == 4;
    res &= s4.all() == true;

    shape s5{0b1111_shape};
    res &= std::ranges::size(s5) == 4;
    res &= s5.all() == true;

    return res;
}

TEST(shape, ctr)
{
    EXPECT_TRUE((std::is_default_constructible_v<shape>));
    EXPECT_TRUE((std::is_nothrow_default_constructible_v<shape>));
    EXPECT_TRUE((std::is_copy_constructible_v<shape>));
    EXPECT_TRUE((std::is_trivially_copy_constructible_v<shape>));
    EXPECT_TRUE((std::is_nothrow_copy_constructible_v<shape>));
    EXPECT_TRUE((std::is_move_constructible_v<shape>));
    EXPECT_TRUE((std::is_trivially_move_constructible_v<shape>));
    EXPECT_TRUE((std::is_nothrow_move_constructible_v<shape>));
    EXPECT_TRUE((std::is_copy_assignable_v<shape>));
    EXPECT_TRUE((std::is_trivially_copy_assignable_v<shape>));
    EXPECT_TRUE((std::is_nothrow_copy_assignable_v<shape>));
    EXPECT_TRUE((std::is_move_assignable_v<shape>));
    EXPECT_TRUE((std::is_trivially_move_assignable_v<shape>));
    EXPECT_TRUE((std::is_nothrow_move_assignable_v<shape>));

    constexpr bool constexpr_ctr = construction_test();
    EXPECT_TRUE(constexpr_ctr);
    bool ctr = construction_test();
    EXPECT_TRUE(ctr);
}

constexpr bool size_test()
{
    bool res{true};
    res = res && (std::ranges::size(shape{ungapped{1}}) == 1u);
    res = res && (std::ranges::size(shape{ungapped{30}}) == 30u);
    res = res && (std::ranges::size(0b11_shape) == 2u);
    res = res && (std::ranges::size(0b10101_shape) == 5u);

    return res;
}

TEST(shape, size)
{
    constexpr bool constexpr_size = size_test();
    EXPECT_TRUE(constexpr_size);
    bool size = size_test();
    EXPECT_TRUE(size);
}
