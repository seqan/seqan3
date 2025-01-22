// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <forward_list>
#include <ranges>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/pretty_printing.hpp>
#include <seqan3/utility/views/interleave.hpp>
#include <seqan3/utility/views/type_reduce.hpp>

using seqan3::operator""_dna4;

TEST(view_interleave, basic)
{
    std::string u{"FOOBARBAXBAT"};
    std::string i{"in"};
    size_t s = 3;
    std::string cmp{"FOOinBARinBAXinBAT"};
    std::string cmp_rev{"TABni"};
    size_t cmpsize = 18;

    // pipe notation
    // explicitly call seqan3::views::type_reduce
    auto v0 = seqan3::views::type_reduce(u) | seqan3::views::interleave(s, seqan3::views::type_reduce(i));
    EXPECT_RANGE_EQ(cmp, v0);
    EXPECT_EQ(cmpsize, v0.size());
    // don't call seqan3::views::type_reduce
    auto v1 = u | seqan3::views::interleave(s, i);
    EXPECT_RANGE_EQ(cmp, v1);

    // function notation
    // explicitly call seqan3::views::type_reduce
    auto v2{seqan3::views::interleave(seqan3::views::type_reduce(u), s, seqan3::views::type_reduce(i))};
    EXPECT_RANGE_EQ(cmp, v2);
    // don't call seqan3::views::type_reduce
    auto v3{seqan3::views::interleave(u, s, i)};
    EXPECT_RANGE_EQ(cmp, v3);

    // combinability
    // explicitly call seqan3::views::type_reduce
    auto v4 = seqan3::views::type_reduce(u) | seqan3::views::interleave(s, seqan3::views::type_reduce(i))
            | std::views::reverse | std::views::take(5);
    EXPECT_RANGE_EQ(cmp_rev, v4);
    // don't call seqan3::views::type_reduce
    auto v5 = u | seqan3::views::interleave(s, i) | std::views::reverse | std::views::take(5);
    EXPECT_RANGE_EQ(cmp_rev, v5);
}

TEST(view_interleave, concepts)
{
    // random_access_range, viewable_range, sized_range
    std::string u{"FOOBARBAXBAT"};
    std::string i{"in"};
    size_t s = 3;
    auto v1 = seqan3::detail::view_interleave(seqan3::views::type_reduce(u), s, seqan3::views::type_reduce(i));

    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE((std::ranges::output_range<decltype(v1), char>));

    EXPECT_FALSE(std::ranges::contiguous_range<decltype(v1)>);

    // forward_range, viewable_range
    std::forward_list<seqan3::dna4> u2{'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4};
    seqan3::dna4_vector i2{'G'_dna4};
    auto v2 = seqan3::views::interleave(u2, s, i2);
    EXPECT_TRUE(std::ranges::input_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::view<decltype(v2)>);

    EXPECT_FALSE(std::ranges::forward_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::contiguous_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v2)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v2), seqan3::dna4>));
}

TEST(view_interleave, chunk_join)
{
    std::forward_list<seqan3::dna4> u{'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4, 'A'_dna4};
    seqan3::dna4_vector i{'G'_dna4};
    size_t s = 2;

    seqan3::dna4_vector cmp{"AAGAAGAA"_dna4};

    auto v1 = seqan3::views::interleave(u, s, i);

    auto it = v1.begin();
    for (auto c : cmp)
    {
        EXPECT_EQ(c, *it);
        it++;
    }
}
