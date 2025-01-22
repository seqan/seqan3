// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream/optional.hpp>
#include <seqan3/core/debug_stream/range.hpp>
#include <seqan3/core/debug_stream/tuple.hpp>
#include <seqan3/core/debug_stream/variant.hpp>
#include <seqan3/test/pretty_printing.hpp>

using seqan3::operator""_dna4;
using namespace std::string_literals;

// Returns a string as gtest would print the given value.
// Output might change due to changes in gtest.
auto gtest_str = [](auto && v)
{
    return ::testing::PrintToString(v);
};

// Returns a string as seqan3 would print the given value.
auto debug_str = [](auto && v)
{
    std::stringstream sstream{};
    seqan3::debug_stream_type dstream{sstream};
    dstream << v;
    return sstream.str();
};

TEST(pretty_printing, char)
{
    EXPECT_EQ(gtest_str('a'), "'a' (97, 0x61)"s);
    EXPECT_EQ(debug_str('a'), "a"s);
}

TEST(pretty_printing, charX_t)
{
    EXPECT_EQ(gtest_str(char8_t{5}), "U+0005"s);
    EXPECT_EQ(gtest_str(char16_t{5}), "U+0005"s);
    EXPECT_EQ(gtest_str(char32_t{5}), "U+0005"s);
    // EXPECT_EQ(debug_str(char8_t{5}), "U+0005"s);
    // EXPECT_EQ(debug_str(char16_t{5}), "U+0005"s);
    // EXPECT_EQ(debug_str(char32_t{5}), "U+0005"s);
}

TEST(pretty_printing, cstring)
{
    EXPECT_EQ(gtest_str("test"), "\"test\""s);
    EXPECT_EQ(debug_str("test"), "test"s);
}

TEST(pretty_printing, tuple)
{
    EXPECT_EQ(gtest_str(std::make_tuple<int, int>(42, -10)), "(42, -10)"s);
    EXPECT_EQ(debug_str(std::make_tuple<int, int>(42, -10)), "(42,-10)"s);
}

TEST(pretty_printing, variant)
{
    EXPECT_EQ(gtest_str(std::variant<int>{0}), "('int(index = 0)' with value 0)"s);
    EXPECT_EQ(debug_str(std::variant<int>{0}), "0"s);
}

TEST(pretty_printing, optional)
{
    EXPECT_EQ(gtest_str(std::optional<int>{}), "(nullopt)"s);
    EXPECT_EQ(debug_str(std::optional<int>{}), "<VALUELESS_OPTIONAL>"s);

    EXPECT_EQ(gtest_str(std::nullopt), "(nullopt)"s);
    EXPECT_EQ(debug_str(std::nullopt), "<VALUELESS_OPTIONAL>"s);
}

TEST(pretty_printing, vector)
{
    EXPECT_EQ(gtest_str(std::vector<std::vector<int>>{{0, 1}, {2, 3}, {1, 2}, {0}}), "[[0,1],[2,3],[1,2],[0]]"s);
    EXPECT_EQ(debug_str(std::vector<std::vector<int>>{{0, 1}, {2, 3}, {1, 2}, {0}}), "[[0,1],[2,3],[1,2],[0]]"s);
}

TEST(pretty_printing, dna)
{
    EXPECT_EQ(gtest_str('G'_dna4), "G"s);
    EXPECT_EQ(debug_str('G'_dna4), "G"s);
}

TEST(pretty_printing, dna_sequence)
{
    EXPECT_EQ(gtest_str("ACGTCGA"_dna4), "ACGTCGA"s);
    EXPECT_EQ(debug_str("ACGTCGA"_dna4), "ACGTCGA"s);

    std::vector dna_2d{"AC"_dna4, "GT"_dna4, "CG"_dna4, "A"_dna4};
    EXPECT_EQ(gtest_str(dna_2d), "[AC,GT,CG,A]"s);
    EXPECT_EQ(debug_str(dna_2d), "[AC,GT,CG,A]"s);

    std::vector dna_3d{std::vector{"AC"_dna4, "GT"_dna4}, std::vector{"CG"_dna4, "A"_dna4}};
    EXPECT_EQ(gtest_str(dna_3d), "[[AC,GT],[CG,A]]"s);
    EXPECT_EQ(debug_str(dna_3d), "[[AC,GT],[CG,A]]"s);
}

TEST(pretty_printing, dna_tuple)
{
    auto dna_tuple = std::make_tuple('A'_dna4, 'C'_dna4);
    EXPECT_EQ(gtest_str(dna_tuple), "(A, C)"s);
    EXPECT_EQ(debug_str(dna_tuple), "(A,C)"s);

    auto dna_sequence_tuple = std::make_tuple("AC"_dna4, "GT"_dna4);
    EXPECT_EQ(gtest_str(dna_sequence_tuple), "(AC, GT)"s);
    EXPECT_EQ(debug_str(dna_sequence_tuple), "(AC,GT)"s);
}

namespace seqan3::detail
{

struct my_type
{
    std::string str;

    friend bool operator==(my_type const & lhs, my_type const & rhs)
    {
        return lhs.str == rhs.str;
    }

    friend bool operator!=(my_type const & lhs, my_type const & rhs)
    {
        return lhs.str != rhs.str;
    }
};

} // namespace seqan3::detail

namespace seqan3
{

template <typename char_t, typename my_type>
    requires std::same_as<std::decay_t<my_type>, detail::my_type>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, my_type && m)
{
    s << m.str;
    return s;
}

} // namespace seqan3

TEST(pretty_printing, seqan3_detail)
{
    // seqan3 types should always produce the same result
    EXPECT_EQ(gtest_str(seqan3::detail::my_type{"HALLO"s}), "HALLO"s);
    EXPECT_EQ(debug_str(seqan3::detail::my_type{"HALLO"s}), "HALLO"s);
}

namespace seqan3
{

struct your_type : public detail::my_type
{};

template <typename char_t, typename your_type>
    requires std::same_as<std::decay_t<your_type>, ::seqan3::your_type>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, your_type && m)
{
    s << m.str;
    return s;
}

} // namespace seqan3

TEST(pretty_printing, seqan3_detail_inherit)
{
    // seqan3 types should always produce the same result
    EXPECT_EQ(gtest_str(seqan3::your_type{"HALLO"s}), "HALLO"s);
    EXPECT_EQ(debug_str(seqan3::your_type{"HALLO"s}), "HALLO"s);
}
