// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>

#include <seqan3/search/kmer_index/shape.hpp>

using seqan3::operator""_shape;

constexpr bool construction_test()
{
    bool res{true};

    seqan3::shape s1{seqan3::bin_literal{0b1011}};
    res &= std::ranges::size(s1) == 4;
    res &= s1.all() == false;

    seqan3::shape s2{0b1011_shape};
    res &= std::ranges::size(s2) == 4;
    res &= s2.all() == false;

    seqan3::shape s3{seqan3::ungapped{3}};
    res &= std::ranges::size(s3) == 3;
    res &= s3.all() == true;

    seqan3::shape s4{seqan3::bin_literal{0b1111}};
    res &= std::ranges::size(s4) == 4;
    res &= s4.all() == true;

    seqan3::shape s5{0b1111_shape};
    res &= std::ranges::size(s5) == 4;
    res &= s5.all() == true;

    return res;
}

TEST(shape, ctr)
{
    EXPECT_TRUE((std::is_default_constructible_v<seqan3::shape>));
    EXPECT_TRUE((std::is_nothrow_default_constructible_v<seqan3::shape>));
    EXPECT_TRUE((std::is_copy_constructible_v<seqan3::shape>));
    EXPECT_TRUE((std::is_trivially_copy_constructible_v<seqan3::shape>));
    EXPECT_TRUE((std::is_nothrow_copy_constructible_v<seqan3::shape>));
    EXPECT_TRUE((std::is_move_constructible_v<seqan3::shape>));
    EXPECT_TRUE((std::is_trivially_move_constructible_v<seqan3::shape>));
    EXPECT_TRUE((std::is_nothrow_move_constructible_v<seqan3::shape>));
    EXPECT_TRUE((std::is_copy_assignable_v<seqan3::shape>));
    EXPECT_TRUE((std::is_trivially_copy_assignable_v<seqan3::shape>));
    EXPECT_TRUE((std::is_nothrow_copy_assignable_v<seqan3::shape>));
    EXPECT_TRUE((std::is_move_assignable_v<seqan3::shape>));
    EXPECT_TRUE((std::is_trivially_move_assignable_v<seqan3::shape>));
    EXPECT_TRUE((std::is_nothrow_move_assignable_v<seqan3::shape>));

    constexpr bool constexpr_ctr = construction_test();
    EXPECT_TRUE(constexpr_ctr);
    bool ctr = construction_test();
    EXPECT_TRUE(ctr);
}

constexpr bool size_test()
{
    bool res{true};
    res = res && (std::ranges::size(seqan3::shape{seqan3::ungapped{1}}) == 1u);
    res = res && (std::ranges::size(seqan3::shape{seqan3::ungapped{30}}) == 30u);
    res = res && (std::ranges::size(0b11_shape) == 2u);
    res = res && (std::ranges::size(0b1'0101_shape) == 5u);

    return res;
}

TEST(shape, size)
{
    constexpr bool constexpr_size = size_test();
    EXPECT_TRUE(constexpr_size);
    bool size = size_test();
    EXPECT_TRUE(size);
}
