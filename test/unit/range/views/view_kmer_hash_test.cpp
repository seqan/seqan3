// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <forward_list>
#include <list>
#include <type_traits>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/range/views/to.hpp>

#include <gtest/gtest.h>

using namespace seqan3;

class kmer_hash_test : public ::testing::Test
{
protected:
    using result_t = std::vector<size_t>;

    static constexpr auto ungapped_view = views::kmer_hash(ungapped{3});
    static constexpr auto gapped_view = views::kmer_hash(0b101_shape);

    std::vector<dna4> text1{"AAAAA"_dna4};
    std::vector<dna4> const ctext1{"AAAAA"_dna4};
    result_t ungapped1{0,0,0};
    result_t gapped1{0,0,0};

    std::vector<dna4> text2{"ACGTAGC"_dna4};
    std::vector<dna4> const ctext2{"ACGTAGC"_dna4};
    result_t ungapped2{6,27,44,50,9};
    result_t gapped2{2, 7, 8, 14, 1};

    std::vector<dna4> text3{"AC"_dna4};
    std::vector<dna4> const ctext3{"AC"_dna4};
    result_t ungapped3{};
    result_t gapped3{ungapped3};

    bitcompressed_vector<dna4> text4{"ACGTAGC"_dna4};
    bitcompressed_vector<dna4> const ctext4{"ACGTAGC"_dna4};
    result_t ungapped4{ungapped2};
    result_t gapped4{gapped2};

    std::list<dna4> text5{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4};
    std::list<dna4> const ctext5{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4};
    result_t ungapped5{ungapped2};
    result_t gapped5{gapped2};

    std::forward_list<dna4> text6{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4};
    std::forward_list<dna4> const ctext6{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'A'_dna4, 'G'_dna4, 'C'_dna4};
    result_t ungapped6{ungapped2};
    result_t gapped6{gapped2};
};


TEST_F(kmer_hash_test, concepts)
{
    auto v1 = text1 | ungapped_view;
    EXPECT_TRUE(std::ranges::input_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::random_access_range<decltype(v1)>);
    EXPECT_TRUE(std::ranges::view<decltype(v1)>);
    EXPECT_TRUE(std::ranges::sized_range<decltype(v1)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v1)>);
    EXPECT_TRUE(const_iterable_range<decltype(v1)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v1), size_t>));

    auto v2 = text5 | ungapped_view;
    EXPECT_TRUE(std::ranges::input_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::bidirectional_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v2)>);
    EXPECT_TRUE(std::ranges::view<decltype(v2)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v2)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v2)>);
    EXPECT_TRUE(const_iterable_range<decltype(v2)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v2), size_t>));

    auto v3 = text6 | ungapped_view;
    EXPECT_TRUE(std::ranges::input_range<decltype(v3)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v3)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v3)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v3)>);
    EXPECT_TRUE(std::ranges::view<decltype(v3)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v3)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v3)>);
    EXPECT_TRUE(const_iterable_range<decltype(v3)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v3), size_t>));
}

TEST_F(kmer_hash_test, ungapped)
{
    EXPECT_EQ(ungapped1, text1 | ungapped_view | views::to<result_t>);
    EXPECT_EQ(ungapped2, text2 | ungapped_view | views::to<result_t>);
    EXPECT_EQ(ungapped3, text3 | ungapped_view | views::to<result_t>);
    EXPECT_EQ(ungapped4, text4 | ungapped_view | views::to<result_t>);
    EXPECT_EQ(ungapped5, text5 | ungapped_view | views::to<result_t>);
    EXPECT_EQ(ungapped6, text6 | ungapped_view | views::to<result_t>);
}

TEST_F(kmer_hash_test, gapped)
{
    EXPECT_EQ(gapped1, text1 | gapped_view | views::to<result_t>);
    EXPECT_EQ(gapped2, text2 | gapped_view | views::to<result_t>);
    EXPECT_EQ(gapped3, text3 | gapped_view | views::to<result_t>);
    EXPECT_EQ(gapped4, text4 | gapped_view | views::to<result_t>);
    EXPECT_EQ(gapped5, text5 | gapped_view | views::to<result_t>);
    EXPECT_EQ(gapped6, text6 | gapped_view | views::to<result_t>);
}

TEST_F(kmer_hash_test, const_ungapped)
{
    EXPECT_EQ(ungapped1, ctext1 | ungapped_view | views::to<result_t>);
    EXPECT_EQ(ungapped2, ctext2 | ungapped_view | views::to<result_t>);
    EXPECT_EQ(ungapped3, ctext3 | ungapped_view | views::to<result_t>);
    EXPECT_EQ(ungapped4, ctext4 | ungapped_view | views::to<result_t>);
    EXPECT_EQ(ungapped5, ctext5 | ungapped_view | views::to<result_t>);
    EXPECT_EQ(ungapped6, ctext6 | ungapped_view | views::to<result_t>);
}

TEST_F(kmer_hash_test, const_gapped)
{
    EXPECT_EQ(gapped1, ctext1 | gapped_view | views::to<result_t>);
    EXPECT_EQ(gapped2, ctext2 | gapped_view | views::to<result_t>);
    EXPECT_EQ(gapped3, ctext3 | gapped_view | views::to<result_t>);
    EXPECT_EQ(gapped4, ctext4 | gapped_view | views::to<result_t>);
    EXPECT_EQ(gapped5, ctext5 | gapped_view | views::to<result_t>);
    EXPECT_EQ(gapped6, ctext6 | gapped_view | views::to<result_t>);
}

TEST_F(kmer_hash_test, combinability)
{
    auto stop_at_t = views::take_until([] (dna4 const x) { return x == 'T'_dna4; });
    EXPECT_EQ(result_t{6}, text2 | stop_at_t | ungapped_view | views::to<result_t>);
    EXPECT_EQ(result_t{6}, text5 | stop_at_t | ungapped_view | views::to<result_t>);
    EXPECT_EQ(result_t{6}, text6 | stop_at_t | ungapped_view | views::to<result_t>);

    EXPECT_EQ(ungapped2 | std::views::reverse | views::to<result_t>,
              text2 | ungapped_view | std::views::reverse | views::to<result_t>);

    EXPECT_EQ(gapped2 | std::views::reverse | views::to<result_t>,
              text2 | gapped_view | std::views::reverse | views::to<result_t>);

    EXPECT_EQ(ungapped5 | std::views::reverse | views::to<result_t>,
              text5 | ungapped_view | std::views::reverse | views::to<result_t>);

    EXPECT_EQ(gapped5 | std::views::reverse | views::to<result_t>,
              text5 | gapped_view | std::views::reverse | views::to<result_t>);
}

TEST_F(kmer_hash_test, invalid_sizes)
{
    EXPECT_THROW(text1 | views::kmer_hash(ungapped{33}), std::invalid_argument);
    EXPECT_THROW(text1 | std::views::reverse | views::kmer_hash(ungapped{33}), std::invalid_argument);
}
