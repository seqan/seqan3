// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <sstream>

#include <seqan3/utility/type_list/detail/type_list_algorithm.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

//-----------------------------------------------------------------------------
// all_of
//-----------------------------------------------------------------------------

struct is_integral_fn
{
    bool operator()(...)
    {
        return false;
    }

    template <typename identity_t>
        requires std::integral<typename identity_t::type>
    bool operator()(identity_t)
    {
        return true;
    }
};

TEST(pack_algorithm, all_of_in_type_list)
{
    EXPECT_TRUE(seqan3::detail::all_of<seqan3::type_list<>>(is_integral_fn{}));
    EXPECT_TRUE((seqan3::detail::all_of<seqan3::type_list<int8_t, int16_t, uint32_t>>(is_integral_fn{})));
    EXPECT_FALSE((seqan3::detail::all_of<seqan3::type_list<int8_t, int16_t, uint32_t, float>>(is_integral_fn{})));
}

//-----------------------------------------------------------------------------
// for_each
//-----------------------------------------------------------------------------

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
