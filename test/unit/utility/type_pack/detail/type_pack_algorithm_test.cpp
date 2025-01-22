// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <sstream>

#include <seqan3/utility/type_pack/detail/type_pack_algorithm.hpp>

//-----------------------------------------------------------------------------
// all_of
//-----------------------------------------------------------------------------

auto is_value_type_integral = [](auto value)
{
    return std::is_integral_v<decltype(value)>;
};

TEST(type_pack_algorithm, all_of)
{
    EXPECT_TRUE(seqan3::detail::all_of(is_value_type_integral));
    EXPECT_TRUE((seqan3::detail::all_of(is_value_type_integral, int8_t{}, int16_t{}, uint32_t{})));
    EXPECT_FALSE((seqan3::detail::all_of(is_value_type_integral, int8_t{}, int16_t{}, uint32_t{}, float{})));
}

//-----------------------------------------------------------------------------
// for_each
//-----------------------------------------------------------------------------

TEST(type_pack_algorithm, for_each)
{
    int i = 0;
    auto fn = [&i](int arg)
    {
        EXPECT_EQ(i, arg);
        ++i;
    };

    seqan3::detail::for_each(fn);
    EXPECT_EQ(i, 0);
    seqan3::detail::for_each(fn, 0);
    EXPECT_EQ(i, 1);
    seqan3::detail::for_each(fn, 1, 2);
    EXPECT_EQ(i, 3);
    seqan3::detail::for_each(fn, 3, 4, 5);
    EXPECT_EQ(i, 6);
}

struct alphabet
{
    char chr;
};

TEST(type_pack_algorithm, for_each2)
{
    std::stringstream stream{};

    auto fn = [&stream](auto const & arg)
    {
        if constexpr (std::is_same_v<decltype(arg), alphabet const &>)
            stream << arg.chr << ";";
        else
            stream << arg << ";";
    };

    seqan3::detail::for_each(fn);
    seqan3::detail::for_each(fn, 0);
    seqan3::detail::for_each(fn, 1.0, '2');
    seqan3::detail::for_each(fn, "3;4", -5, alphabet{'C'});

    EXPECT_EQ(stream.str(), "0;1;2;3;4;-5;C;");
}
