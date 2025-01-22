// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <ranges>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/alphabet/views/trim_quality.hpp>
#include <seqan3/test/expect_range_eq.hpp>

using seqan3::operator""_dna5;
using seqan3::operator""_phred42;

TEST(view_trim, standalone)
{
    using namespace std::literals;

    std::vector<seqan3::phred42> vec{"II?5+"_phred42};

    // trim by phred_value
    EXPECT_RANGE_EQ("II?5"_phred42, vec | seqan3::views::trim_quality(20u));

    // trim by quality character
    EXPECT_RANGE_EQ("II"_phred42, vec | seqan3::views::trim_quality('I'_phred42));

    // function syntax
    EXPECT_RANGE_EQ("II?5"_phred42, seqan3::views::trim_quality(vec, '5'_phred42));

    // combinability
    EXPECT_RANGE_EQ("II?5"sv, seqan3::views::trim_quality(vec, 20u) | seqan3::views::to_char);
}

TEST(view_trim, qualified)
{
    using namespace std::literals;

    std::vector<seqan3::dna5q> vec{{'A'_dna5, 'I'_phred42},
                                   {'G'_dna5, 'I'_phred42},
                                   {'G'_dna5, '?'_phred42},
                                   {'A'_dna5, '5'_phred42},
                                   {'T'_dna5, '+'_phred42}};
    std::vector<seqan3::dna5q> cmp1{{'A'_dna5, 'I'_phred42},
                                    {'G'_dna5, 'I'_phred42},
                                    {'G'_dna5, '?'_phred42},
                                    {'A'_dna5, '5'_phred42}};
    std::vector<seqan3::dna5q> cmp2{{'A'_dna5, 'I'_phred42}, {'G'_dna5, 'I'_phred42}};

    // trim by phred_value
    EXPECT_RANGE_EQ(cmp1, vec | seqan3::views::trim_quality(20u));

    // trim by quality character
    EXPECT_RANGE_EQ(cmp2, vec | seqan3::views::trim_quality(seqan3::dna5q{'C'_dna5, 'I'_phred42}));

    // function syntax
    EXPECT_RANGE_EQ(cmp1, seqan3::views::trim_quality(vec, 20u));

    // combinability
    EXPECT_RANGE_EQ("AGGA"sv, seqan3::views::trim_quality(vec, 20u) | seqan3::views::to_char);
}

TEST(view_trim, concepts)
{
    std::vector<seqan3::dna5q> vec{};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), seqan3::dna5q>));
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);

    auto v1 = vec | seqan3::views::trim_quality(20u);
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v1), seqan3::dna5q>));
    EXPECT_TRUE(!std::ranges::sized_range<decltype(v1)>);
}
