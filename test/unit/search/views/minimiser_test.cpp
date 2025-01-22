// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <forward_list>
#include <list>
#include <type_traits>

#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/io/views/detail/take_until_view.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "../../range/iterator_test_template.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_shape;
using result_t = std::vector<size_t>;

static inline constexpr auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{4});
static inline constexpr auto rev_kmer_view =
    seqan3::views::complement | std::views::reverse | kmer_view | std::views::reverse;
static inline constexpr auto gapped_kmer_view = seqan3::views::kmer_hash(0b1001_shape);
static inline constexpr auto rev_gapped_kmer_view =
    seqan3::views::complement | std::views::reverse | seqan3::views::kmer_hash(0b1001_shape) | std::views::reverse;
static inline constexpr auto minimiser_view1 = seqan3::views::minimiser(1); // kmer_size == window_size
static inline constexpr auto minimiser_no_rev_view = seqan3::views::minimiser(5);

using iterator_type =
    std::ranges::iterator_t<decltype(std::declval<seqan3::dna4_vector &>() | kmer_view | minimiser_no_rev_view)>;
using two_ranges_iterator_type = std::ranges::iterator_t<decltype(seqan3::detail::minimiser_view{
    std::declval<seqan3::dna4_vector &>() | kmer_view,
    std::declval<seqan3::dna4_vector &>() | rev_kmer_view,
    5})>;

template <>
struct iterator_fixture<iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4})) vec = text | kmer_view;
    result_t expected_range{26, 97, 27};

    decltype(seqan3::views::minimiser(seqan3::views::kmer_hash(text, seqan3::ungapped{4}), 5)) test_range =
        seqan3::views::minimiser(vec, 5);
};

template <>
struct iterator_fixture<two_ranges_iterator_type> : public ::testing::Test
{
    using iterator_tag = std::forward_iterator_tag;
    static constexpr bool const_iterable = true;

    seqan3::dna4_vector text{"ACGGCGACGTTTAG"_dna4};
    using kmer_hash_view_t = decltype(seqan3::views::kmer_hash(text, seqan3::ungapped{4}));

    kmer_hash_view_t vec = kmer_view(text);
    result_t expected_range{26, 97, 27, 6, 1};

    using reverse_kmer_hash_view_t = decltype(rev_kmer_view(text));

    using test_range_t = decltype(seqan3::detail::minimiser_view{std::declval<kmer_hash_view_t &>(),
                                                                 std::declval<reverse_kmer_hash_view_t &>(),
                                                                 5});
    test_range_t test_range = seqan3::detail::minimiser_view{vec, rev_kmer_view(text), 5};
};

using test_types = ::testing::Types<iterator_type, two_ranges_iterator_type>;
INSTANTIATE_TYPED_TEST_SUITE_P(iterator_fixture, iterator_fixture, test_types, );

template <typename T>
class minimiser_view_properties_test : public ::testing::Test
{};

using underlying_range_types = ::testing::Types<std::vector<seqan3::dna4>,
                                                std::vector<seqan3::dna4> const,
                                                seqan3::bitpacked_sequence<seqan3::dna4>,
                                                seqan3::bitpacked_sequence<seqan3::dna4> const,
                                                std::list<seqan3::dna4>,
                                                std::list<seqan3::dna4> const,
                                                std::forward_list<seqan3::dna4>,
                                                std::forward_list<seqan3::dna4> const>;
TYPED_TEST_SUITE(minimiser_view_properties_test, underlying_range_types, );

class minimiser_test : public ::testing::Test
{
protected:
    std::vector<seqan3::dna4> text1{"AAAAAAAAAAAAAAAAAAA"_dna4};
    std::vector<seqan3::dna4> text1_short{"AAAAAA"_dna4};
    result_t result1{0, 0, 0}; // Same result for ungapped and gapped
    result_t result1_short{0}; // windows_size == text_size, same result for ungapped and gapped

    std::vector<seqan3::dna4> too_short_text{"AC"_dna4};

    std::vector<seqan3::dna4> text3{"ACGGCGACGTTTAG"_dna4};
    result_t result3_ungapped{26, 97, 27, 6, 1};  // ACGG, CGAC, ACGT, aacg, aaac - lowercase for reverse complement
    result_t result3_gapped{2, 5, 3, 2, 1};       // A--G, C--C, A--T, a--g, a--c - "-" for gap
    result_t result3_ungapped_no_rev{26, 97, 27}; // ACGG, CGAC, ACGT
    result_t result3_gapped_no_rev{2, 5, 3};      // A--G, C--C-, A--T "-" for gap
    result_t result3_ungapped_stop{26, 97};       // For stop at first T
    result_t result3_gapped_stop{2, 5};           // For stop at first T
    result_t result3_start{1};                    // For start at second A, ungapped and gapped the same
    result_t result3_ungapped_no_rev_start{27};   // For start at second A
    result_t result3_gapped_no_rev_start{3};      // For start at second A
};

template <typename adaptor_t>
void compare_types(adaptor_t v)
{
    EXPECT_TRUE(std::ranges::input_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::forward_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::bidirectional_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::random_access_range<decltype(v)>);
    EXPECT_TRUE(std::ranges::view<decltype(v)>);
    EXPECT_FALSE(std::ranges::sized_range<decltype(v)>);
    EXPECT_FALSE(std::ranges::common_range<decltype(v)>);
    EXPECT_TRUE(seqan3::const_iterable_range<decltype(v)>);
    EXPECT_FALSE((std::ranges::output_range<decltype(v), size_t>));
}

TYPED_TEST(minimiser_view_properties_test, concepts)
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
                   'G'_dna4}; // ACGTCGACGTTTAG

    auto v = text | kmer_view | minimiser_no_rev_view;
    compare_types(v);
    auto v2 = seqan3::detail::minimiser_view{text | kmer_view, text | kmer_view, 5};

    if constexpr (std::ranges::bidirectional_range<TypeParam>) // excludes forward_list
    {
        auto v3 = seqan3::detail::minimiser_view{text | kmer_view, text | rev_kmer_view, 5};
        compare_types(v3);
    }
}

TYPED_TEST(minimiser_view_properties_test, different_inputs_kmer_hash)
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
                   'G'_dna4};             // ACGTCGACGTTTAG
    result_t ungapped{27, 97, 27, 6, 1};  // ACGT, CGAC, ACGT, aacg, aaac - lowercase for reverse comp.
    result_t gapped{3, 5, 3, 2, 1};       // A--T, C--C, A--T, a--g, a--c - "-" for gap
    result_t ungapped_no_rev{27, 97, 27}; // ACGT, CGAC, ACGT
    result_t gapped_no_rev{3, 5, 3};      // A--T, C--C, A--T - "-" for gap
    EXPECT_RANGE_EQ(ungapped_no_rev, text | kmer_view | minimiser_no_rev_view);
    EXPECT_RANGE_EQ(gapped_no_rev, text | gapped_kmer_view | minimiser_no_rev_view);

    if constexpr (std::ranges::bidirectional_range<TypeParam>) // excludes forward_list
    {
        EXPECT_RANGE_EQ(ungapped, (seqan3::detail::minimiser_view{text | kmer_view, text | rev_kmer_view, 5}));
        EXPECT_RANGE_EQ(gapped,
                        (seqan3::detail::minimiser_view{text | gapped_kmer_view, text | rev_gapped_kmer_view, 5}));
    }
}

TEST_F(minimiser_test, ungapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1, (seqan3::detail::minimiser_view{text1 | kmer_view, text1 | rev_kmer_view, 5}));
    EXPECT_RANGE_EQ(result1, text1 | kmer_view | minimiser_no_rev_view);
    EXPECT_THROW(text1_short | kmer_view | minimiser_view1, std::invalid_argument);
    auto empty_view = seqan3::detail::minimiser_view{too_short_text | kmer_view, too_short_text | rev_kmer_view, 5};
    EXPECT_TRUE(std::ranges::empty(empty_view));
    auto empty_view2 = too_short_text | kmer_view | minimiser_no_rev_view;
    EXPECT_TRUE(std::ranges::empty(empty_view2));
    EXPECT_RANGE_EQ(result3_ungapped, (seqan3::detail::minimiser_view{text3 | kmer_view, text3 | rev_kmer_view, 5}));
    EXPECT_RANGE_EQ(result3_ungapped_no_rev, text3 | kmer_view | minimiser_no_rev_view);
}

TEST_F(minimiser_test, gapped_kmer_hash)
{
    EXPECT_RANGE_EQ(result1,
                    (seqan3::detail::minimiser_view{text1 | gapped_kmer_view, text1 | rev_gapped_kmer_view, 5}));
    EXPECT_RANGE_EQ(result1, text1 | gapped_kmer_view | minimiser_no_rev_view);
    EXPECT_THROW(text1_short | gapped_kmer_view | minimiser_view1, std::invalid_argument);
    auto empty_view =
        seqan3::detail::minimiser_view{too_short_text | gapped_kmer_view, too_short_text | rev_gapped_kmer_view, 5};
    EXPECT_TRUE(std::ranges::empty(empty_view));
    auto empty_view2 = too_short_text | gapped_kmer_view | minimiser_no_rev_view;
    EXPECT_TRUE(std::ranges::empty(empty_view2));
    EXPECT_RANGE_EQ(result3_gapped,
                    (seqan3::detail::minimiser_view{text3 | gapped_kmer_view, text3 | rev_gapped_kmer_view, 5}));
    EXPECT_RANGE_EQ(result3_gapped_no_rev, text3 | gapped_kmer_view | minimiser_no_rev_view);
}

TEST_F(minimiser_test, window_too_big)
{
    EXPECT_RANGE_EQ(result1_short, text1 | kmer_view | seqan3::views::minimiser(20));
    EXPECT_RANGE_EQ(result1_short, text1 | gapped_kmer_view | seqan3::views::minimiser(20));
    EXPECT_RANGE_EQ(result1_short, (seqan3::detail::minimiser_view{text1 | kmer_view, text1 | rev_kmer_view, 20}));
    EXPECT_RANGE_EQ(result1_short,
                    (seqan3::detail::minimiser_view{text1 | gapped_kmer_view, text1 | rev_gapped_kmer_view, 20}));
}

TEST_F(minimiser_test, combinability)
{
    auto stop_at_t = std::views::take_while(
        [](seqan3::dna4 const x)
        {
            return x != 'T'_dna4;
        });
    EXPECT_RANGE_EQ(result3_ungapped_stop, text3 | stop_at_t | kmer_view | minimiser_no_rev_view);
    EXPECT_RANGE_EQ(result3_gapped_stop, text3 | stop_at_t | gapped_kmer_view | minimiser_no_rev_view);

    EXPECT_RANGE_EQ(
        result3_ungapped_stop,
        (seqan3::detail::minimiser_view{text3 | stop_at_t | kmer_view, text3 | stop_at_t | rev_kmer_view, 5}));
    EXPECT_RANGE_EQ(result3_gapped_stop,
                    (seqan3::detail::minimiser_view{text3 | stop_at_t | gapped_kmer_view,
                                                    text3 | stop_at_t | rev_gapped_kmer_view,
                                                    5}));

    auto start_at_a = std::views::drop(6);
    EXPECT_RANGE_EQ(
        result3_start,
        (seqan3::detail::minimiser_view{text3 | start_at_a | kmer_view, text3 | start_at_a | rev_kmer_view, 5}));
    EXPECT_RANGE_EQ(result3_start,
                    (seqan3::detail::minimiser_view{text3 | start_at_a | gapped_kmer_view,
                                                    text3 | start_at_a | rev_gapped_kmer_view,
                                                    5}));
}

TEST_F(minimiser_test, non_arithmetic_value)
{
    // just compute the minimizer directly on the alphabet
    EXPECT_RANGE_EQ("ACACA"_dna4, text3 | minimiser_no_rev_view);
}

TEST_F(minimiser_test, two_ranges_unequal_size)
{
    EXPECT_THROW((seqan3::detail::minimiser_view{text1 | kmer_view, text3 | rev_kmer_view, 5}), std::invalid_argument);
}

TEST_F(minimiser_test, iterator_base)
{
    std::vector<size_t> hashes{3, 6, 5, 4, 8, 4, 4, 2, 5, 4};
    using hash_iterator_t = typename std::vector<size_t>::iterator;

    hash_iterator_t const hash_begin = hashes.begin();
    hash_iterator_t const hash_end = hashes.end();

    size_t const window_size = 5;
    auto const minimiser = seqan3::views::minimiser(hashes, window_size);
    auto minimiser_it = minimiser.begin();

    hash_iterator_t const hash_first_window_end = minimiser_it.base();

    // hash iterator is at position window_size - 1 and
    // points to the last element in the window
    // index:   0, 1, 2, 3, 4,  5, 6, 7, 8, 9
    // hashes: [3, 6, 5, 4, 8], 4, 4, 2, 5, 4
    //                      ^
    EXPECT_EQ(minimiser_it.base() - hash_first_window_end, 0u);    // window start position
    EXPECT_EQ(minimiser_it.base() - hash_begin, window_size - 1u); // window end position

    // After incrementing, points to the last element in the new window
    // index:  0,  1, 2, 3, 4, 5,  6, 7, 8, 9
    // hashes: 3, [6, 5, 4, 8, 4], 4, 2, 5, 4
    //                         ^
    ++minimiser_it;
    EXPECT_EQ(minimiser_it.base() - hash_first_window_end, 1u); // window start position
    EXPECT_EQ(minimiser_it.base() - hash_begin, 5u);            // window end position

    // After incrementing, points to the last element in the new window
    // index:  0, 1, 2,  3, 4, 5, 6, 7,  8, 9
    // hashes: 3, 6, 5, [4, 8, 4, 4, 2], 5, 4
    //                               ^
    ++minimiser_it;
    EXPECT_EQ(minimiser_it.base() - hash_first_window_end, 3u); // window start position
    EXPECT_EQ(minimiser_it.base() - hash_begin, 7u);            // window end position

    // if minimiser reached end, underlying iterator reached end too
    ++minimiser_it;
    EXPECT_EQ(minimiser_it, minimiser.end());
    EXPECT_EQ(minimiser_it.base(), hash_end);
}
