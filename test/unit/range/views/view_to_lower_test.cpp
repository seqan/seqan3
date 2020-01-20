// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <iostream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/to_lower.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

TEST(view_to_lower, basic)
{
    std::string input_string {"IAmADnaString"};
    std::string cmp {"iamadnastring"};

    // pipe notation string
    std::string s(input_string | views::to_lower | views::to<std::string>);
    EXPECT_EQ(cmp, s);

    // custom conversion operator
    std::string s2(views::to_lower(input_string) | views::to<std::string>);
    EXPECT_EQ(cmp, s2);
}

TEST(view_to_lower, combinability)
{
    std::string input_string {"IAmADnaString"};
    std::string cmp{"gnirtsandamai"};

    std::vector<dna5> dna_vec {"AGGCGT"_dna5};
    std::string cmp2{"aggcgt"};

   // output combinability
    std::string s(input_string | views::to_lower | std::views::reverse | views::to<std::string>);
    EXPECT_EQ(cmp, s);

    // input combinability
    std::string s2(dna_vec | views::to_char | views::to_lower | views::to<std::string>);
    EXPECT_EQ(cmp2, s2);
}

TEST(view_to_lower, deep)
{
    std::vector<std::string> input_vec{"IAmADnaString", "IAmAProteinString"};
    std::vector<std::string> cmp{"iamadnastring", "iamaproteinstring"};

    std::vector<std::string> s(input_vec | views::to_lower | views::to<std::vector<std::string>>);
    EXPECT_EQ(cmp, s);
}

TEST(view_to_lower, concepts)
{
    std::string input_string{"AEIOU"};
    std::string & input_string_ref = input_string;
    auto lower_view = input_string | views::to_lower;

    // Required
    EXPECT_TRUE(std::ranges::input_range<decltype(input_string)>);
    EXPECT_TRUE(std::ranges::viewable_range<decltype(input_string_ref)>);

    // Preserved
    EXPECT_EQ(std::ranges::input_range<decltype(input_string)>,
              std::ranges::input_range<decltype(lower_view)>);
    EXPECT_EQ(std::ranges::forward_range<decltype(input_string)>,
              std::ranges::forward_range<decltype(lower_view)>);
    EXPECT_EQ(std::ranges::bidirectional_range<decltype(input_string)>,
              std::ranges::bidirectional_range<decltype(lower_view)>);
    EXPECT_EQ(std::ranges::random_access_range<decltype(input_string)>,
              std::ranges::random_access_range<decltype(lower_view)>);
    EXPECT_EQ(std::ranges::random_access_range<decltype(input_string)>,
              std::ranges::random_access_range<decltype(lower_view)>);
    EXPECT_EQ(std::ranges::viewable_range<decltype(input_string_ref)>,
              std::ranges::viewable_range<decltype(lower_view)>);
    EXPECT_EQ(std::ranges::sized_range<decltype(input_string)>,
              std::ranges::sized_range<decltype(lower_view)>);
    EXPECT_EQ(std::ranges::common_range<decltype(input_string)>,
              std::ranges::common_range<decltype(lower_view)>);
    EXPECT_EQ(const_iterable_range<decltype(input_string)>,
              const_iterable_range<decltype(lower_view)>);
    EXPECT_TRUE((std::same_as<std::remove_reference_t<reference_t<decltype(input_string)>>,
                           std::remove_reference_t<reference_t<decltype(lower_view)>>>));

    // Guaranteed
    EXPECT_TRUE(std::ranges::viewable_range<decltype(lower_view)>);
    EXPECT_TRUE(std::ranges::view<decltype(lower_view)>);

    // Lost
    EXPECT_FALSE((std::ranges::output_range<decltype(lower_view), char>));
    EXPECT_FALSE(std::ranges::contiguous_range<decltype(lower_view)>);
}
