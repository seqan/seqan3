// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>
#include <ranges>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/mask/masked.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/range/concept.hpp>
#include <seqan3/utility/views/elements.hpp>
#include <seqan3/utility/views/zip.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_phred42;

TEST(view_get, basic)
{
    // TODO remove const-ness from input vector once alphabet_proxy's complement doesnt cause ICE
    std::vector<seqan3::dna4q> const qv{{'A'_dna4, '!'_phred42},
                                        {'C'_dna4, '"'_phred42},
                                        {'G'_dna4, '#'_phred42},
                                        {'T'_dna4, '$'_phred42}};

    //functor
    EXPECT_RANGE_EQ("ACGT"_dna4, seqan3::views::elements<0>(qv));
    EXPECT_RANGE_EQ("!\"#$"_phred42, seqan3::views::elements<1>(qv));

    // pipe notation
    EXPECT_RANGE_EQ("ACGT"_dna4, qv | seqan3::views::elements<0>);
    EXPECT_RANGE_EQ("!\"#$"_phred42, qv | seqan3::views::elements<1>);

    // combinability
    EXPECT_RANGE_EQ("TGCA"_dna4, qv | seqan3::views::elements<0> | seqan3::views::complement);
}

TEST(view_get, advanced)
{
    // TODO remove const-ness from input vector once alphabet_proxy inherits it's alphabet
    std::vector<seqan3::qualified<seqan3::masked<seqan3::dna4>, seqan3::phred42>> const t{
        {{'A'_dna4, seqan3::mask::masked}, '!'_phred42},
        {{'C'_dna4, seqan3::mask::unmasked}, '"'_phred42},
        {{'G'_dna4, seqan3::mask::masked}, '#'_phred42},
        {{'T'_dna4, seqan3::mask::unmasked}, '$'_phred42}};

    // functor notation
    std::vector<seqan3::masked<seqan3::dna4>> expected_sequence{{'A'_dna4, seqan3::mask::masked},
                                                                {'C'_dna4, seqan3::mask::unmasked},
                                                                {'G'_dna4, seqan3::mask::masked},
                                                                {'T'_dna4, seqan3::mask::unmasked}};

    EXPECT_RANGE_EQ(expected_sequence, seqan3::views::elements<0>(t));
    EXPECT_RANGE_EQ("!\"#$"_phred42, seqan3::views::elements<1>(t));

    EXPECT_RANGE_EQ("ACGT"_dna4, seqan3::views::elements<0>(seqan3::views::elements<0>(t)));

    // pipe notation
    EXPECT_RANGE_EQ(expected_sequence, t | seqan3::views::elements<0>);
    EXPECT_RANGE_EQ("!\"#$"_phred42, t | seqan3::views::elements<1>);

    EXPECT_RANGE_EQ("ACGT"_dna4, t | seqan3::views::elements<0> | seqan3::views::elements<0>);

    // combinability
    EXPECT_RANGE_EQ(expected_sequence | std::views::reverse, t | seqan3::views::elements<0> | std::views::reverse);
    EXPECT_RANGE_EQ("$#\"!"_phred42, t | seqan3::views::elements<1> | std::views::reverse);

    EXPECT_RANGE_EQ("TGCA"_dna4, t | seqan3::views::elements<0> | seqan3::views::elements<0> | std::views::reverse);
}

TEST(view_get, pair_range)
{
    std::vector<std::pair<int, int>> pair_range{{0, 1}, {1, 2}, {2, 3}, {3, 4}};

    // functor notation
    EXPECT_RANGE_EQ((std::vector{0, 1, 2, 3}), seqan3::views::elements<0>(pair_range));
    EXPECT_RANGE_EQ((std::vector{1, 2, 3, 4}), seqan3::views::elements<1>(pair_range));

    // pipe notation
    EXPECT_RANGE_EQ((std::vector{0, 1, 2, 3}), pair_range | seqan3::views::elements<0>);
    EXPECT_RANGE_EQ((std::vector{1, 2, 3, 4}), pair_range | seqan3::views::elements<1>);
}

TEST(view_get, tuple_range)
{
    std::vector<std::tuple<int, int>> tuple_range{{0, 1}, {1, 2}, {2, 3}, {3, 4}};

    // functor notation
    EXPECT_RANGE_EQ((std::vector{0, 1, 2, 3}), seqan3::views::elements<0>(tuple_range));
    EXPECT_RANGE_EQ((std::vector{1, 2, 3, 4}), seqan3::views::elements<1>(tuple_range));

    // pipe notation
    EXPECT_RANGE_EQ((std::vector{0, 1, 2, 3}), tuple_range | seqan3::views::elements<0>);
    EXPECT_RANGE_EQ((std::vector{1, 2, 3, 4}), tuple_range | seqan3::views::elements<1>);
}

TEST(view_get, concepts)
{
    std::vector<std::tuple<int, int>> vec{{0, 1}, {0, 1}, {0, 1}, {0, 1}, {0, 1}};
    EXPECT_TRUE(std::ranges::input_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(vec)>);
    EXPECT_FALSE(std::ranges::view<decltype(vec)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(vec)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(vec)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(vec)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(vec), std::tuple<int, int>>));

    auto v1 = vec | seqan3::views::elements<0>;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), std::tuple<int, int>>));
    EXPECT_TRUE((std::ranges::output_range<decltype(v1), int>));
}

// https://github.com/seqan/seqan3/issues/745
TEST(view_get, nested_zip_view)
{
    std::vector vec1{0, 1, 2, 3, 4, 5, 6, 7, 8, 9};

    auto get_view = seqan3::views::zip(seqan3::views::zip(vec1, vec1), vec1) | seqan3::views::elements<0>;

    for (auto && elem : get_view)
        std::get<0>(elem) = -1;

    EXPECT_EQ(vec1[0], -1);
}
