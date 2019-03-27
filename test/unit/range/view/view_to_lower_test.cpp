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
#include <seqan3/range/view/to_lower.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

TEST(view_to_lower, basic)
{
    std::string input_string {"IAmADnaString"};
    std::string cmp {"iamadnastring"};

    // pipe notation string
    std::string s(input_string | view::to_lower);
    EXPECT_EQ(cmp, s);

    // custom conversion operator
    std::string s2(view::to_lower(input_string));
    EXPECT_EQ(cmp, s2);
}

TEST(view_to_lower, combinability)
{
    std::string input_string {"IAmADnaString"};
    std::string cmp{"gnirtsandamai"};

    std::vector<dna5> dna_vec {"AGGCGT"_dna5};
    std::string cmp2{"aggcgt"};

   // output combinability
    std::string s(input_string | view::to_lower | std::view::reverse);
    EXPECT_EQ(cmp, s);

    // input combinability
    std::string s2(dna_vec | view::to_char | view::to_lower);
    EXPECT_EQ(cmp2, s2);
}

TEST(view_to_lower, deep)
{
    std::vector<std::string> input_vec{"IAmADnaString", "IAmAProteinString"};
    std::vector<std::string> cmp{"iamadnastring", "iamaproteinstring"};

    std::vector<std::string> s(input_vec | view::to_lower);
    EXPECT_EQ(cmp, s);
}

TEST(view_to_lower, concepts)
{
    std::string input_string{"AEIOU"};
    std::string & input_string_ref = input_string;
    auto lower_view = input_string | view::to_lower;

    // Required
    EXPECT_TRUE(std::ranges::InputRange<decltype(input_string)>);
    EXPECT_TRUE(std::ranges::ViewableRange<decltype(input_string_ref)>);

    // Preserved
    EXPECT_EQ(std::ranges::InputRange<decltype(input_string)>,
              std::ranges::InputRange<decltype(lower_view)>);
    EXPECT_EQ(std::ranges::ForwardRange<decltype(input_string)>,
              std::ranges::ForwardRange<decltype(lower_view)>);
    EXPECT_EQ(std::ranges::BidirectionalRange<decltype(input_string)>,
              std::ranges::BidirectionalRange<decltype(lower_view)>);
    EXPECT_EQ(std::ranges::RandomAccessRange<decltype(input_string)>,
              std::ranges::RandomAccessRange<decltype(lower_view)>);
    EXPECT_EQ(std::ranges::RandomAccessRange<decltype(input_string)>,
              std::ranges::RandomAccessRange<decltype(lower_view)>);
    EXPECT_EQ(std::ranges::ViewableRange<decltype(input_string_ref)>,
              std::ranges::ViewableRange<decltype(lower_view)>);
    EXPECT_EQ(std::ranges::SizedRange<decltype(input_string)>,
              std::ranges::SizedRange<decltype(lower_view)>);
    EXPECT_EQ(std::ranges::CommonRange<decltype(input_string)>,
              std::ranges::CommonRange<decltype(lower_view)>);
    EXPECT_EQ(const_iterable_concept<decltype(input_string)>,
              const_iterable_concept<decltype(lower_view)>);
    EXPECT_TRUE((std::Same<std::remove_reference_t<reference_t<decltype(input_string)>>,
                           std::remove_reference_t<reference_t<decltype(lower_view)>>>));

    // Guaranteed
    EXPECT_TRUE(std::ranges::ViewableRange<decltype(lower_view)>);
    EXPECT_TRUE(std::ranges::View<decltype(lower_view)>);

    // Lost
    EXPECT_FALSE((std::ranges::OutputRange<decltype(lower_view), char>));
    EXPECT_FALSE(std::ranges::ContiguousRange<decltype(lower_view)>);
}
