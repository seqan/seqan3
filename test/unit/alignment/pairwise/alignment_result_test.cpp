// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_rna5;

using aligned_seq_type = std::vector<seqan3::gapped<seqan3::dna4>>;

template <typename T>
struct alignment_result_test : public ::testing::Test
{};

using alignment_result_test_types =
    ::testing::Types<seqan3::detail::alignment_result_value_type<uint32_t,
                                                                 uint32_t,
                                                                 int32_t,
                                                                 std::pair<size_t, size_t>,
                                                                 std::pair<size_t, size_t>,
                                                                 std::pair<aligned_seq_type, aligned_seq_type>>,
                     seqan3::detail::alignment_result_value_type<uint32_t,
                                                                 uint32_t,
                                                                 int32_t,
                                                                 std::pair<size_t, size_t>,
                                                                 std::pair<size_t, size_t>,
                                                                 std::tuple<aligned_seq_type, aligned_seq_type>>,
                     seqan3::detail::alignment_result_value_type<uint32_t,
                                                                 uint32_t,
                                                                 int32_t,
                                                                 std::pair<size_t, size_t>,
                                                                 std::pair<size_t, size_t>,
                                                                 std::vector<aligned_seq_type>>,
                     seqan3::detail::alignment_result_value_type<uint32_t,
                                                                 uint32_t,
                                                                 float,
                                                                 std::pair<size_t, size_t>,
                                                                 std::pair<size_t, size_t>,
                                                                 std::pair<aligned_seq_type, aligned_seq_type>>,
                     seqan3::detail::alignment_result_value_type<uint32_t,
                                                                 uint32_t,
                                                                 float,
                                                                 std::pair<size_t, size_t>,
                                                                 std::pair<size_t, size_t>,
                                                                 std::tuple<aligned_seq_type, aligned_seq_type>>,
                     seqan3::detail::alignment_result_value_type<uint32_t,
                                                                 uint32_t,
                                                                 float,
                                                                 std::pair<size_t, size_t>,
                                                                 std::pair<size_t, size_t>,
                                                                 std::vector<aligned_seq_type>>>;

TYPED_TEST_SUITE(alignment_result_test, alignment_result_test_types, );

TYPED_TEST(alignment_result_test, type_specialisation)
{
    EXPECT_TRUE((seqan3::detail::is_type_specialisation_of_v<TypeParam, seqan3::detail::alignment_result_value_type>));
}

TYPED_TEST(alignment_result_test, constructor)
{
    EXPECT_TRUE((std::is_default_constructible_v<seqan3::alignment_result<TypeParam>>));
    EXPECT_TRUE((std::is_copy_constructible_v<seqan3::alignment_result<TypeParam>>));
    EXPECT_TRUE((std::is_move_constructible_v<seqan3::alignment_result<TypeParam>>));
    EXPECT_TRUE((std::is_copy_assignable_v<seqan3::alignment_result<TypeParam>>));
    EXPECT_TRUE((std::is_move_assignable_v<seqan3::alignment_result<TypeParam>>));
    EXPECT_TRUE((std::is_destructible_v<seqan3::alignment_result<TypeParam>>));
}

TYPED_TEST(alignment_result_test, get_id)
{
    using id1_t = decltype(std::declval<TypeParam>().sequence1_id);
    using id2_t = decltype(std::declval<TypeParam>().sequence2_id);
    aligned_seq_type seq{'A'_dna4, 'T'_dna4, seqan3::gap{}, 'C'_dna4, seqan3::gap{}, seqan3::gap{}, 'A'_dna4};

    {
        seqan3::alignment_result<TypeParam> tmp{TypeParam{1u, 2u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence1_id()), id1_t>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence2_id()), id2_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).sequence1_id()), id1_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).sequence2_id()), id2_t>));
        EXPECT_EQ(tmp.sequence1_id(), 1u);
        EXPECT_EQ(tmp.sequence2_id(), 2u);
    }

    {
        seqan3::alignment_result<TypeParam> const tmp{TypeParam{1u, 2u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.sequence1_id(), 1u);
        EXPECT_EQ(tmp.sequence2_id(), 2u);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence1_id()), id1_t>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence2_id()), id2_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).sequence1_id()), id1_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).sequence2_id()), id2_t>));
    }
}

TYPED_TEST(alignment_result_test, get_score)
{
    using score_t = decltype(std::declval<TypeParam>().score);
    aligned_seq_type seq{'A'_dna4, 'T'_dna4, seqan3::gap{}, 'C'_dna4, seqan3::gap{}, seqan3::gap{}, 'A'_dna4};

    {
        seqan3::alignment_result<TypeParam> tmp{TypeParam{1u, 2u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.score(), 0);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.score()), score_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).score()), score_t>));
    }

    {
        seqan3::alignment_result<TypeParam> const tmp{TypeParam{1u, 2u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.score(), 0);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.score()), score_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).score()), score_t>));
    }
}

TYPED_TEST(alignment_result_test, end_positions)
{
    aligned_seq_type seq{'A'_dna4, 'T'_dna4, seqan3::gap{}, 'C'_dna4, seqan3::gap{}, seqan3::gap{}, 'A'_dna4};

    {
        seqan3::alignment_result<TypeParam> tmp{TypeParam{1u, 2u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.sequence1_end_position(), 10ul);
        EXPECT_EQ(tmp.sequence2_end_position(), 10ul);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence1_end_position()), size_t>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence2_end_position()), size_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).sequence1_end_position()), size_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).sequence2_end_position()), size_t>));
    }

    {
        seqan3::alignment_result<TypeParam> const tmp{TypeParam{1u, 2u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.sequence1_end_position(), 10ul);
        EXPECT_EQ(tmp.sequence2_end_position(), 10ul);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence1_end_position()), size_t>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence2_end_position()), size_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).sequence1_end_position()), size_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).sequence2_end_position()), size_t>));
    }
}

TYPED_TEST(alignment_result_test, begin_positions)
{
    aligned_seq_type seq{'A'_dna4, 'T'_dna4, seqan3::gap{}, 'C'_dna4, seqan3::gap{}, seqan3::gap{}, 'A'_dna4};

    {
        seqan3::alignment_result<TypeParam> tmp{TypeParam{1u, 2u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.sequence1_begin_position(), 0ul);
        EXPECT_EQ(tmp.sequence2_begin_position(), 0ul);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence1_begin_position()), size_t>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence2_begin_position()), size_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).sequence1_begin_position()), size_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).sequence2_begin_position()), size_t>));
    }

    {
        seqan3::alignment_result<TypeParam> const tmp{TypeParam{1u, 2u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.sequence1_begin_position(), 0ul);
        EXPECT_EQ(tmp.sequence2_begin_position(), 0ul);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence1_begin_position()), size_t>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence2_begin_position()), size_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).sequence1_begin_position()), size_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).sequence2_begin_position()), size_t>));
    }
}

TYPED_TEST(alignment_result_test, alignment)
{
    using namespace std::literals;

    using alignment_t = decltype(std::declval<TypeParam>().alignment);
    aligned_seq_type seq{'A'_dna4, 'T'_dna4, seqan3::gap{}, 'C'_dna4, seqan3::gap{}, seqan3::gap{}, 'A'_dna4};

    {
        seqan3::alignment_result<TypeParam> tmp{TypeParam{1u, 2u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.alignment(), (alignment_t{seq, seq}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.alignment()), alignment_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).alignment()), alignment_t const &>));
    }

    {
        seqan3::alignment_result<TypeParam> const tmp{TypeParam{1u, 2u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.alignment(), (alignment_t{seq, seq}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.alignment()), alignment_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).alignment()), alignment_t const &>));
    }

    if constexpr (seqan3::tuple_like<alignment_t>)
    {
        seqan3::alignment_result<TypeParam> tmp{TypeParam{1u, 2u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_RANGE_EQ(std::get<0>(tmp.alignment()) | std::views::all | seqan3::views::to_char, "AT-C--A"sv);
        EXPECT_RANGE_EQ(std::get<1>(tmp.alignment()) | std::views::all | seqan3::views::to_char, "AT-C--A"sv);
    }
    else
    {
        seqan3::alignment_result<TypeParam> tmp{TypeParam{1u, 2u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_RANGE_EQ(tmp.alignment()[0] | std::views::all | seqan3::views::to_char, "AT-C--A"sv);
        EXPECT_RANGE_EQ(tmp.alignment()[1] | std::views::all | seqan3::views::to_char, "AT-C--A"sv);
    }
}

TEST(alignment_result_test, reduced_type)
{
    {
        seqan3::detail::alignment_result_value_type tr{2u, 4u, 5};
        seqan3::alignment_result tmp(tr);
        EXPECT_EQ(tmp.sequence1_id(), 2u);
        EXPECT_EQ(tmp.sequence2_id(), 4u);
        EXPECT_EQ(tmp.score(), 5);
    }

    {
        seqan3::detail::alignment_result_value_type tr{2, 4u, 5.0f, std::pair<int, int>{1, -1}};
        seqan3::alignment_result tmp(tr);
        EXPECT_EQ(tmp.sequence1_id(), 2);
        EXPECT_EQ(tmp.sequence2_id(), 4u);
        EXPECT_FLOAT_EQ(tmp.score(), 5.0f);
        EXPECT_EQ((tmp.sequence1_end_position()), 1);
        EXPECT_EQ((tmp.sequence2_end_position()), -1);
    }

    {
        seqan3::detail::alignment_result_value_type tr{2,
                                                       4u,
                                                       5.0f,
                                                       std::pair<int, int>{1, -1},
                                                       std::pair<int, int>{10, -10}};
        seqan3::alignment_result tmp(tr);
        EXPECT_EQ(tmp.sequence1_id(), 2);
        EXPECT_EQ(tmp.sequence2_id(), 4u);
        EXPECT_FLOAT_EQ(tmp.score(), 5.0f);
        EXPECT_EQ((tmp.sequence1_end_position()), 1);
        EXPECT_EQ((tmp.sequence2_end_position()), -1);
        EXPECT_EQ((tmp.sequence1_begin_position()), 10);
        EXPECT_EQ((tmp.sequence2_begin_position()), -10);
    }
}

TEST(alignment_result_test, type_deduction)
{
    {
        using coord_t = std::pair<int, int>;
        std::vector<seqan3::gapped<seqan3::rna5>>
            seq{'A'_rna5, 'U'_rna5, seqan3::gap{}, 'C'_rna5, seqan3::gap{}, seqan3::gap{}, 'A'_rna5};

        seqan3::detail::alignment_result_value_type tr{2, 4, 5.0, coord_t{1, -1}, coord_t{10, -10}, seq};
        seqan3::alignment_result tmp(tr);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence1_id()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence2_id()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.score()), double>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence1_end_position()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence2_end_position()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence1_begin_position()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence2_begin_position()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.alignment()), std::vector<seqan3::gapped<seqan3::rna5>> const &>));

        EXPECT_EQ(tmp.sequence1_id(), 2);
        EXPECT_EQ(tmp.sequence2_id(), 4);
        EXPECT_DOUBLE_EQ(tmp.score(), 5.0);
        EXPECT_EQ((tmp.sequence1_end_position()), 1);
        EXPECT_EQ((tmp.sequence2_end_position()), -1);
        EXPECT_EQ((tmp.sequence1_begin_position()), 10);
        EXPECT_EQ((tmp.sequence2_begin_position()), -10);
        EXPECT_EQ(tmp.alignment(), seq);
    }

    {
        using coord_t = std::pair<int, int>;

        seqan3::detail::alignment_result_value_type tr{2, 4, 5.0, coord_t{1, -1}};
        seqan3::alignment_result tmp(tr);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence1_id()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence2_id()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.score()), double>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence1_end_position()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence2_end_position()), int>));
    }

    {
        seqan3::detail::alignment_result_value_type tr{2, 3, 5.0};
        seqan3::alignment_result tmp(tr);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence1_id()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.sequence2_id()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.score()), double>));
    }
}

TEST(alignment_result_test, empty_type)
{
    seqan3::detail::alignment_result_value_type tr{};
    seqan3::alignment_result tmp(tr);
    // members must not be accessed
}

TEST(alignment_result_test, access_result_value_type)
{
    seqan3::detail::alignment_result_value_type result_value{2u, 4u, 5};
    seqan3::alignment_result result(result_value);

    using result_value_t = typename seqan3::detail::alignment_result_value_type_accessor<decltype(result)>::type;
    EXPECT_SAME_TYPE(result_value_t, decltype(result_value));
}
