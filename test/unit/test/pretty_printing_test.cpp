// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>
#include <optional>
#include <ostream>
#include <variant>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/pretty_printing.hpp>

using seqan3::operator""_dna4;
using namespace std::string_literals;

// Returns a string as gtest would print the given value.
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

TEST(pretty_printing, gtest_output)
{
    // warning: these outputs can change due to changes in gtest
    EXPECT_EQ(gtest_str('a'), "'a' (97, 0x61)"s);
    EXPECT_EQ(debug_str('a'), "a"s);

    EXPECT_EQ(gtest_str(std::make_tuple<int, int>(42, -10)), "(42, -10)"s);
    EXPECT_EQ(debug_str(std::make_tuple<int, int>(42, -10)), "(42,-10)"s);
}

TEST(pretty_printing, std_output)
{
    EXPECT_EQ(gtest_str(std::vector<std::vector<int>>{{0,1}, {2,3}, {1,2}, {0}}), "[[0,1],[2,3],[1,2],[0]]"s);
    EXPECT_EQ(debug_str(std::vector<std::vector<int>>{{0,1}, {2,3}, {1,2}, {0}}), "[[0,1],[2,3],[1,2],[0]]"s);

    EXPECT_EQ(gtest_str(std::variant<int>{0}), "0"s);
    EXPECT_EQ(debug_str(std::variant<int>{0}), "0"s);

    EXPECT_EQ(gtest_str(std::optional<int>{}), "<VALUELESS_OPTIONAL>"s);
    EXPECT_EQ(debug_str(std::optional<int>{}), "<VALUELESS_OPTIONAL>"s);

    EXPECT_EQ(gtest_str(std::nullopt), "<VALUELESS_OPTIONAL>"s);
    EXPECT_EQ(debug_str(std::nullopt), "<VALUELESS_OPTIONAL>"s);
}

TEST(pretty_printing, seqan3_output)
{
    // seqan3 types should always produce the same result
    EXPECT_EQ(gtest_str('G'_dna4), "G"s);
    EXPECT_EQ(debug_str('G'_dna4), "G"s);
    EXPECT_EQ('G'_dna4, 'G'_dna4); // change value to test

    EXPECT_EQ(gtest_str("ACGTCGA"_dna4), "ACGTCGA"s);
    EXPECT_EQ(debug_str("ACGTCGA"_dna4), "ACGTCGA"s);
    EXPECT_EQ("ACGTCGA"_dna4, "ACGTCGA"_dna4); // change value to test

    std::vector dna_set1{"AC"_dna4, "GT"_dna4, "CG"_dna4, "A"_dna4};
    std::vector dna_set2{"AC"_dna4, "GT"_dna4, "CG"_dna4, "A"_dna4};
    EXPECT_EQ(gtest_str(dna_set1), "[AC,GT,CG,A]"s);
    EXPECT_EQ(debug_str(dna_set1), "[AC,GT,CG,A]"s);
    EXPECT_EQ(dna_set1, dna_set2); // change value to test

    std::vector dna_sets1{std::vector{"AC"_dna4, "GT"_dna4}, std::vector{"CG"_dna4, "A"_dna4}};
    std::vector dna_sets2{std::vector{"AC"_dna4, "GT"_dna4}, std::vector{"CG"_dna4, "A"_dna4}};
    EXPECT_EQ(gtest_str(dna_sets1), "[[AC,GT],[CG,A]]"s);
    EXPECT_EQ(debug_str(dna_sets1), "[[AC,GT],[CG,A]]"s);
    EXPECT_EQ(dna_sets1, dna_sets2); // change value to test
}

TEST(pretty_printing, gtest_mixed_seqan3_output)
{
    // warning: these outputs can change due to changes in gtest
    auto dna_tuple1 = std::make_tuple('A'_dna4, 'C'_dna4);
    auto dna_tuple2 = std::make_tuple('A'_dna4, 'C'_dna4);
    EXPECT_EQ(gtest_str(dna_tuple1), "(A, C)"s);
    EXPECT_EQ(debug_str(dna_tuple1), "(A,C)"s);
    EXPECT_EQ(dna_tuple1, dna_tuple2); // change value to test

    auto dna_sequence_tuple1 = std::make_tuple("AC"_dna4, "GT"_dna4);
    auto dna_sequence_tuple2 = std::make_tuple("AC"_dna4, "GT"_dna4);
    EXPECT_EQ(gtest_str(dna_sequence_tuple1), "(AC, GT)"s);
    EXPECT_EQ(debug_str(dna_sequence_tuple1), "(AC,GT)"s);
    EXPECT_EQ(dna_sequence_tuple1, dna_sequence_tuple2); // change value to test
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

template <typename my_type, typename char_t>
    requires std::same_as<std::decay_t<my_type>, detail::my_type>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, my_type && m)
{
    s << m.str;
    return s;
}

} // namespace seqan3

TEST(pretty_printing, seqan3_detail_output)
{
    // seqan3 types should always produce the same result
    EXPECT_EQ(gtest_str(seqan3::detail::my_type{"HALLO"s}), "HALLO"s);
    EXPECT_EQ(debug_str(seqan3::detail::my_type{"HALLO"s}), "HALLO"s);
    EXPECT_EQ(seqan3::detail::my_type{"HALLO"}, seqan3::detail::my_type{"HALLO"}); // change value to test
}

namespace seqan3
{

struct your_type : public detail::my_type
{};

template <typename your_type, typename char_t>
    requires std::same_as<std::decay_t<your_type>, ::seqan3::your_type>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, your_type && m)
{
    s << m.str;
    return s;
}

} // namespace seqan3

TEST(pretty_printing, seqan3_detail_mixed_seqan3_output)
{
    // seqan3 types should always produce the same result
    EXPECT_EQ(gtest_str(seqan3::your_type{"HALLO"s}), "HALLO"s);
    EXPECT_EQ(debug_str(seqan3::your_type{"HALLO"s}), "HALLO"s);
    EXPECT_EQ(seqan3::your_type{"HALLO"}, seqan3::your_type{"HALLO"}); // change value to test
}
