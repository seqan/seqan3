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
#include <seqan3/range/views/to_upper.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/range/concept.hpp>

#ifdef SEQAN3_DEPRECATED_310
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

using seqan3::operator""_dna5;

TEST(view_to_upper, basic)
{
    using namespace std::literals;

    std::string input_string {"IAmADnaString"};

    // pipe notation string
    EXPECT_RANGE_EQ("IAMADNASTRING"sv, input_string | seqan3::views::to_upper);

    // custom conversion operator
    EXPECT_RANGE_EQ("IAMADNASTRING"sv, seqan3::views::to_upper(input_string));
}

TEST(view_to_upper, combinability)
{
    using namespace std::literals;

    // output combinability
    std::string input_string {"IAmADnaString"};
    EXPECT_RANGE_EQ("GNIRTSANDAMAI"sv, input_string | seqan3::views::to_upper | std::views::reverse);

    // input combinability
    std::vector<seqan3::dna5> dna_vec {"aggcgt"_dna5};
    EXPECT_RANGE_EQ("AGGCGT"sv, dna_vec | seqan3::views::to_char | seqan3::views::to_upper);
}

TEST(view_to_upper, deep)
{
    using namespace std::literals;

    std::vector<std::string> input_vec{"IAmADnaString", "IAmAProteinString"};

    auto view = input_vec | seqan3::views::to_upper;
    EXPECT_RANGE_EQ("IAMADNASTRING"sv, view[0]);
    EXPECT_RANGE_EQ("IAMAPROTEINSTRING"sv, view[1]);
}

TEST(view_to_upper, concepts)
{
    std::string input_string{"aeiou"};
    std::string & input_string_ref = input_string;
    auto upper_view = input_string | seqan3::views::to_upper;

    // Required
    EXPECT_TRUE(std::ranges::input_range<decltype(input_string)>);
    EXPECT_TRUE(std::ranges::viewable_range<decltype(input_string_ref)>);

    // Preserved
    EXPECT_EQ(std::ranges::input_range<decltype(input_string)>,
              std::ranges::input_range<decltype(upper_view)>);
    EXPECT_EQ(std::ranges::forward_range<decltype(input_string)>,
              std::ranges::forward_range<decltype(upper_view)>);
    EXPECT_EQ(std::ranges::bidirectional_range<decltype(input_string)>,
              std::ranges::bidirectional_range<decltype(upper_view)>);
    EXPECT_EQ(std::ranges::random_access_range<decltype(input_string)>,
              std::ranges::random_access_range<decltype(upper_view)>);
    EXPECT_EQ(std::ranges::random_access_range<decltype(input_string)>,
              std::ranges::random_access_range<decltype(upper_view)>);
    EXPECT_EQ(std::ranges::viewable_range<decltype(input_string_ref)>,
              std::ranges::viewable_range<decltype(upper_view)>);
    EXPECT_EQ(std::ranges::sized_range<decltype(input_string)>,
              std::ranges::sized_range<decltype(upper_view)>);
    EXPECT_EQ(std::ranges::common_range<decltype(input_string)>,
              std::ranges::common_range<decltype(upper_view)>);
    EXPECT_EQ(seqan3::const_iterable_range<decltype(input_string)>,
              seqan3::const_iterable_range<decltype(upper_view)>);
    EXPECT_TRUE((std::same_as<std::remove_reference_t<std::ranges::range_reference_t<decltype(input_string)>>,
                              std::remove_reference_t<std::ranges::range_reference_t<decltype(upper_view)>>>));

    // Guaranteed
    EXPECT_TRUE(std::ranges::viewable_range<decltype(upper_view)>);
    EXPECT_TRUE(std::ranges::view<decltype(upper_view)>);

    // Lost
    EXPECT_FALSE((std::ranges::output_range<decltype(upper_view), char>));
    EXPECT_FALSE(std::ranges::contiguous_range<decltype(upper_view)>);
}
#pragma GCC diagnostic pop
#endif // SEQAN3_DEPRECATED_310
