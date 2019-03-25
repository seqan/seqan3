// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <iostream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/to_upper.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

TEST(view_to_upper, basic)
{
    std::string input_string {"IAmADnaString"};
    std::string cmp {"IAMADNASTRING"};

    // pipe notation string
    std::string s(input_string | view::to_upper);
    EXPECT_EQ(cmp, s);

    // custom conversion operator
    std::string s2(view::to_upper(input_string));
    EXPECT_EQ(cmp, s2);
}

TEST(view_to_upper, combinability)
{
    std::string input_string {"IAmADnaString"};
    std::string cmp{"GNIRTSANDAMAI"};

    std::vector<dna5> dna_vec {"aggcgt"_dna5};
    std::string cmp2{"AGGCGT"};

   // output combinability
    std::string s(input_string | view::to_upper | std::view::reverse);
    EXPECT_EQ(cmp, s);

    // input combinability
    std::string s2(dna_vec | view::to_char | view::to_upper);
    EXPECT_EQ(cmp2, s2);
}

TEST(view_to_upper, deep)
{
    std::vector<std::string> input_vec{"IAmADnaString", "IAmAProteinString"};
    std::vector<std::string> cmp{"IAMADNASTRING", "IAMAPROTEINSTRING"};

    std::vector<std::string> s(input_vec | view::to_upper);
    EXPECT_EQ(cmp, s);
}

TEST(view_to_upper, concepts)
{
    std::string input_string{"aeiou"};
    std::string & input_string_ref = input_string;
    auto upper_view = input_string | view::to_upper;

    // Required
    EXPECT_TRUE(std::ranges::InputRange<decltype(input_string)>);
    EXPECT_TRUE(std::ranges::ViewableRange<decltype(input_string_ref)>);

    // Preserved
    EXPECT_EQ(std::ranges::InputRange<decltype(input_string)>,
              std::ranges::InputRange<decltype(upper_view)>);
    EXPECT_EQ(std::ranges::ForwardRange<decltype(input_string)>,
              std::ranges::ForwardRange<decltype(upper_view)>);
    EXPECT_EQ(std::ranges::BidirectionalRange<decltype(input_string)>,
              std::ranges::BidirectionalRange<decltype(upper_view)>);
    EXPECT_EQ(std::ranges::RandomAccessRange<decltype(input_string)>,
              std::ranges::RandomAccessRange<decltype(upper_view)>);
    EXPECT_EQ(std::ranges::RandomAccessRange<decltype(input_string)>,
              std::ranges::RandomAccessRange<decltype(upper_view)>);
    EXPECT_EQ(std::ranges::ViewableRange<decltype(input_string_ref)>,
              std::ranges::ViewableRange<decltype(upper_view)>);
    EXPECT_EQ(std::ranges::SizedRange<decltype(input_string)>,
              std::ranges::SizedRange<decltype(upper_view)>);
    EXPECT_EQ(std::ranges::CommonRange<decltype(input_string)>,
              std::ranges::CommonRange<decltype(upper_view)>);
    EXPECT_EQ(const_iterable_concept<decltype(input_string)>,
              const_iterable_concept<decltype(upper_view)>);
    EXPECT_TRUE((std::Same<std::remove_reference_t<reference_t<decltype(input_string)>>,
                           std::remove_reference_t<reference_t<decltype(upper_view)>>>));

    // Guaranteed
    EXPECT_TRUE(std::ranges::ViewableRange<decltype(upper_view)>);
    EXPECT_TRUE(std::ranges::View<decltype(upper_view)>);

    // Lost
    EXPECT_FALSE((std::ranges::OutputRange<decltype(upper_view), char>));
    EXPECT_FALSE(std::ranges::ContiguousRange<decltype(upper_view)>);
}
