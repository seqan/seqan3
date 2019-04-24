// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>
#include <ostream>

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

    EXPECT_EQ(gtest_str(std::vector<std::vector<int>>{{0,1}, {2,3}, {1,2}, {0}}),
              "{ { 0, 1 }, { 2, 3 }, { 1, 2 }, { 0 } }"s);
    EXPECT_EQ(debug_str(std::vector<std::vector<int>>{{0,1}, {2,3}, {1,2}, {0}}),
              "[[0,1],[2,3],[1,2],[0]]"s);
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

namespace seqan3::detail
{

struct my_type
{
    std::string str;

    bool operator==(my_type const & other) const
    {
        return str == other.str;
    }

    bool operator!=(my_type const & other) const
    {
        return str != other.str;
    }
};

} // namespace seqan3::detail

namespace seqan3
{

template <typename my_type>
    requires std::Same<std::decay_t<my_type>, detail::my_type>
inline debug_stream_type & operator<<(debug_stream_type & s, my_type && m)
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
