// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/ranges>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/range/views/to_lower.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/range/concept.hpp>

#ifdef SEQAN3_DEPRECATED_310
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

using seqan3::operator""_dna5;

TEST(view_to_lower, basic)
{
    using namespace std::literals;

    std::string input_string {"IAmADnaString"};

    // pipe notation string
    EXPECT_RANGE_EQ("iamadnastring"sv, input_string | seqan3::views::to_lower);

    // custom conversion operator
    EXPECT_RANGE_EQ("iamadnastring"sv, seqan3::views::to_lower(input_string));
}

TEST(view_to_lower, combinability)
{
    using namespace std::literals;

    // output combinability
    std::string input_string {"IAmADnaString"};
    EXPECT_RANGE_EQ("gnirtsandamai"sv, input_string | seqan3::views::to_lower | std::views::reverse);

    // input combinability
    std::vector<seqan3::dna5> dna_vec {"AGGCGT"_dna5};
    EXPECT_RANGE_EQ("aggcgt"sv, dna_vec | seqan3::views::to_char | seqan3::views::to_lower);
}

TEST(view_to_lower, deep)
{
    using namespace std::literals;

    std::vector<std::string> input_vec{"IAmADnaString", "IAmAProteinString"};

    auto view = input_vec | seqan3::views::to_lower;
    EXPECT_RANGE_EQ("iamadnastring"sv, view[0]);
    EXPECT_RANGE_EQ("iamaproteinstring"sv, view[1]);
}

TEST(view_to_lower, concepts)
{
    std::string input_string{"AEIOU"};
    std::string & input_string_ref = input_string;
    auto lower_view = input_string | seqan3::views::to_lower;

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
    EXPECT_EQ(seqan3::const_iterable_range<decltype(input_string)>,
              seqan3::const_iterable_range<decltype(lower_view)>);
    EXPECT_TRUE((std::same_as<std::remove_reference_t<std::ranges::range_reference_t<decltype(input_string)>>,
                              std::remove_reference_t<std::ranges::range_reference_t<decltype(lower_view)>>>));

    // Guaranteed
    EXPECT_TRUE(std::ranges::viewable_range<decltype(lower_view)>);
    EXPECT_TRUE(std::ranges::view<decltype(lower_view)>);

    // Lost
    EXPECT_FALSE((std::ranges::output_range<decltype(lower_view), char>));
    EXPECT_FALSE(std::ranges::contiguous_range<decltype(lower_view)>);
}
#pragma GCC diagnostic pop
#endif // SEQAN3_DEPRECATED_310
