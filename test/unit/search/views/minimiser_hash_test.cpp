// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <forward_list>
#include <list>

#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "../../range/iterator_test_template.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;
using iterator_type = std::ranges::iterator_t<
    decltype(std::declval<seqan3::dna4_vector &>()
             | seqan3::views::minimiser_hash(seqan3::ungapped{4}, seqan3::window_size{8}, seqan3::seed{0}))>;

static constexpr seqan3::shape ungapped_shape = seqan3::ungapped{4};
static constexpr seqan3::shape gapped_shape = 0b1001_shape;
static constexpr auto ungapped_view =
    seqan3::views::minimiser_hash(ungapped_shape, seqan3::window_size{8}, seqan3::seed{0});
static constexpr auto gapped_view =
    seqan3::views::minimiser_hash(gapped_shape, seqan3::window_size{8}, seqan3::seed{0});

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = false;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    result_t expected_range{26, 97, 27, 6, 1};

    using test_range_t = decltype(text | ungapped_view);
    test_range_t test_range = text | ungapped_view;
};

using test_type = ::testing::Types<iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_type, );

template <typename T>
class minimiser_hash_properties_test : public ::testing::Test
{};

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitpacked_sequence<seqan3::dna4>,
                                                seqan3::bitpacked_sequence<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const>;

TYPED_TEST_SUITE(minimiser_hash_properties_test, underlying_range_types, );
class minimiser_hash_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAAAAAAAAAAAAAAA"_dna4};
    std::vector<seqan3::dna4> text1_short{"AAAAAA"_dna4};
    result_t result1{0, 0, 0}; // Same for ungapped and gapped
    result_t ungapped_default_seed{0x8F'3F'73'B5'CF'1C'9A'21, 0x8F'3F'73'B5'CF'1C'9A'21, 0x8F'3F'73'B5'CF'1C'9A'21};
    result_t gapped_default_seed{0x8F'3F'73'B5'CF'1C'9A'D1, 0x8F'3F'73'B5'CF'1C'9A'D1, 0x8F'3F'73'B5'CF'1C'9A'D1};

    std::vector<seqan3::dna4> text2{"AC"_dna4};
    result_t result2{};

    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
    result_t ungapped3{26, 97, 27, 6, 1}; // ACGG, CGAC, ACGT, aacg, aaac
    result_t ungapped_stop_at_t3{26, 97}; // ACGG, CGAC
    result_t gapped3{2, 5, 3, 2, 1};      // A--G, C--C, A--T, a--g, a--c "-" for gap
    result_t gapped_stop_at_t3{2, 5};     // A--G, C--C "-" for gap
};

TYPED_TEST(minimiser_hash_properties_test, different_input_ranges)
{
    TypeParam text{'A'_dna4,
                   'C'_dna4,
                   'G'_dna4,
                   'T'_dna4,
                   'C'_dna4,
                   'G'_dna4,
                   'A'_dna4,
                   'C'_dna4,
                   'G'_dna4,
                   'T'_dna4,
                   'T'_dna4,
                   'T'_dna4,
                   'A'_dna4,
                   'G'_dna4};            // ACGTCGACGTTTAG
    result_t ungapped{27, 97, 27, 6, 1}; // ACGT, CGAC, ACGT, aacg, aaac
    result_t gapped{3, 5, 3, 2, 1};      // A--T, C--C, A--T, a--g, a--c - "-" for gap
    EXPECT_RANGE_EQ(ungapped, text | ungapped_view);
    EXPECT_RANGE_EQ(gapped, text | gapped_view);
}

TEST_F(minimiser_hash_test, ungapped)
{
    EXPECT_RANGE_EQ(result1, text1 | ungapped_view);
    EXPECT_RANGE_EQ(result2, text2 | ungapped_view);
    EXPECT_RANGE_EQ(ungapped3, text3 | ungapped_view);

    auto stop_at_t = std::views::take_while(
        [](seqan3::dna4 const x)
        {
            return x != 'T'_dna4;
        });
    EXPECT_RANGE_EQ(ungapped_stop_at_t3, text3 | stop_at_t | ungapped_view);
}

TEST_F(minimiser_hash_test, gapped)
{
    EXPECT_RANGE_EQ(result1, text1 | gapped_view);
    EXPECT_RANGE_EQ(result2, text2 | gapped_view);
    EXPECT_RANGE_EQ(gapped3, text3 | gapped_view);

    auto stop_at_t = std::views::take_while(
        [](seqan3::dna4 const x)
        {
            return x != 'T'_dna4;
        });
    EXPECT_RANGE_EQ(gapped_stop_at_t3, text3 | stop_at_t | gapped_view);
}

TEST_F(minimiser_hash_test, seed)
{
    EXPECT_RANGE_EQ(ungapped_default_seed,
                    text1 | seqan3::views::minimiser_hash(ungapped_shape, seqan3::window_size{8}));
    EXPECT_RANGE_EQ(gapped_default_seed, text1 | seqan3::views::minimiser_hash(gapped_shape, seqan3::window_size{8}));
}

TEST_F(minimiser_hash_test, shape_bigger_than_window)
{
    EXPECT_THROW(text1 | seqan3::views::minimiser_hash(ungapped_shape, seqan3::window_size{3}), std::invalid_argument);
    EXPECT_THROW(text1 | seqan3::views::minimiser_hash(gapped_shape, seqan3::window_size{3}), std::invalid_argument);
}
