// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <sstream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/detail/pack_algorithm.hpp>
#include <seqan3/core/type_list/type_list.hpp>

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
    EXPECT_TRUE(seqan3::detail::all_of<seqan3::type_list<>>(is_integral_fn{}));
    EXPECT_TRUE((seqan3::detail::all_of<seqan3::type_list<int8_t, int16_t, uint32_t>>(is_integral_fn{})));
    EXPECT_FALSE((seqan3::detail::all_of<seqan3::type_list<int8_t, int16_t, uint32_t, float>>(is_integral_fn{})));
}

TEST(pack_algorithm, all_of_values)
{
    EXPECT_TRUE(seqan3::detail::all_of(is_value_type_integral));
    EXPECT_TRUE((seqan3::detail::all_of(is_value_type_integral, int8_t{}, int16_t{}, uint32_t{})));
    EXPECT_FALSE((seqan3::detail::all_of(is_value_type_integral, int8_t{}, int16_t{}, uint32_t{}, float{})));
}

//-----------------------------------------------------------------------------
// for_each
//-----------------------------------------------------------------------------

TEST(pack_algorithm, for_each_value)
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

TEST(pack_algorithm, for_each_value2)
{
    using seqan3::operator""_dna4;

    std::stringstream stream{};

    auto fn = [&stream](auto const & arg)
    {
        if constexpr(seqan3::alphabet<decltype(arg)>)
            stream << seqan3::to_char(arg) << ";";
        else
            stream << arg << ";";
    };

    seqan3::detail::for_each(fn);
    seqan3::detail::for_each(fn, 0);
    seqan3::detail::for_each(fn, 1.0, '2');
    seqan3::detail::for_each(fn, "3;4", -5, 'C'_dna4);

    EXPECT_EQ(stream.str(), "0;1;2;3;4;-5;C;");
}

template <typename type>
void print_to_stream(std::stringstream & stream, std::type_identity<type>)
{
    if constexpr (std::is_same_v<type, bool>)
        stream << type{0} << ";";
    if constexpr (std::is_same_v<type, uint8_t>)
        stream << unsigned{1} << ";";
    if constexpr (std::is_same_v<type, int8_t>)
        stream << int{-1} << ";";
    if constexpr (std::is_same_v<type, uint16_t>)
        stream << type{2} << ";";
    if constexpr (std::is_same_v<type, int16_t>)
        stream << type{-2} << ";";
    if constexpr (std::is_same_v<type, uint32_t>)
        stream << type{3} << ";";
    if constexpr (std::is_same_v<type, int32_t>)
        stream << type{-3} << ";";
    if constexpr (std::is_same_v<type, uint64_t>)
        stream << type{4} << ";";
    if constexpr (std::is_same_v<type, int64_t>)
        stream << type{-4} << ";";
}

TEST(pack_algorithm, for_each_type_in_type_list)
{
    std::stringstream stream{};

    auto fn = [&stream](auto id)
    {
        print_to_stream(stream, id);
    };

    using types = seqan3::type_list<bool, uint8_t, int8_t, uint16_t, int16_t, uint32_t, int32_t, uint64_t, int64_t>;
    seqan3::detail::for_each<types>(fn);

    EXPECT_EQ(stream.str(), "0;1;-1;2;-2;3;-3;4;-4;");

    stream.str("");
    seqan3::detail::for_each<types>(fn);

    EXPECT_EQ(stream.str(), "0;1;-1;2;-2;3;-3;4;-4;");
}
