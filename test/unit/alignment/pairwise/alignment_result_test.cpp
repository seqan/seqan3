// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <meta/meta.hpp>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/range/views/to.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_rna5;

using aligned_seq_type = std::vector<seqan3::gapped<seqan3::dna4>>;

template <typename T>
struct alignment_result_test : public ::testing::Test
{};

using alignment_result_test_types = ::testing::Types
      <seqan3::detail::alignment_result_value_type<uint32_t, int32_t,
                                                   std::pair<size_t, size_t>, std::pair<size_t, size_t>,
                                                   std::pair<aligned_seq_type, aligned_seq_type>>,
       seqan3::detail::alignment_result_value_type<uint32_t, int32_t,
                                                   std::pair<size_t, size_t>, std::pair<size_t, size_t>,
                                                   std::tuple<aligned_seq_type, aligned_seq_type>>,
       seqan3::detail::alignment_result_value_type<uint32_t, int32_t,
                                                   std::pair<size_t, size_t>, std::pair<size_t, size_t>,
                                                   std::vector<aligned_seq_type>>,
       seqan3::detail::alignment_result_value_type<uint32_t, float,
                                                   std::pair<size_t, size_t>, std::pair<size_t, size_t>,
                                                   std::pair<aligned_seq_type, aligned_seq_type>>,
       seqan3::detail::alignment_result_value_type<uint32_t, float,
                                                   std::pair<size_t, size_t>, std::pair<size_t, size_t>,
                                                   std::tuple<aligned_seq_type, aligned_seq_type>>,
       seqan3::detail::alignment_result_value_type<uint32_t, float,
                                                   std::pair<size_t, size_t>, std::pair<size_t, size_t>,
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
    using id_t = decltype(std::declval<TypeParam>().id);
    aligned_seq_type seq{'A'_dna4, 'T'_dna4, seqan3::gap{}, 'C'_dna4, seqan3::gap{}, seqan3::gap{}, 'A'_dna4};

    {
        seqan3::alignment_result<TypeParam> tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_TRUE((std::is_same_v<decltype(tmp.id()), id_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).id()), id_t>));
        EXPECT_EQ(tmp.id(), 1u);
    }

    {
        seqan3::alignment_result<TypeParam> const tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.id(), 1u);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.id()), id_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).id()), id_t>));
    }
}

TYPED_TEST(alignment_result_test, get_score)
{
    using score_t = decltype(std::declval<TypeParam>().score);
    aligned_seq_type seq{'A'_dna4, 'T'_dna4, seqan3::gap{}, 'C'_dna4, seqan3::gap{}, seqan3::gap{}, 'A'_dna4};

    {
        seqan3::alignment_result<TypeParam> tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.score(), 0);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.score()), score_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).score()), score_t>));
    }

    {
        seqan3::alignment_result<TypeParam> const tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.score(), 0);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.score()), score_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).score()), score_t>));
    }

}

TYPED_TEST(alignment_result_test, back_coordinate)
{
    using coord_t = decltype(std::declval<TypeParam>().back_coordinate);
    aligned_seq_type seq{'A'_dna4, 'T'_dna4, seqan3::gap{}, 'C'_dna4, seqan3::gap{}, seqan3::gap{}, 'A'_dna4};

    {
        seqan3::alignment_result<TypeParam> tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.back_coordinate(), (coord_t{10ul, 10ul}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.back_coordinate()), coord_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).back_coordinate()), coord_t const &>));
    }

    {
        seqan3::alignment_result<TypeParam> const tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.back_coordinate(), (coord_t{10ul, 10ul}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.back_coordinate()), coord_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).back_coordinate()), coord_t const &>));
    }
}

TYPED_TEST(alignment_result_test, front_coordinate)
{
    using coord_t = decltype(std::declval<TypeParam>().front_coordinate);
    aligned_seq_type seq{'A'_dna4, 'T'_dna4, seqan3::gap{}, 'C'_dna4, seqan3::gap{}, seqan3::gap{}, 'A'_dna4};

    {
        seqan3::alignment_result<TypeParam> tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.front_coordinate(), (coord_t{0ul, 0ul}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.front_coordinate()), coord_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).front_coordinate()), coord_t const &>));
    }

    {
        seqan3::alignment_result<TypeParam> const tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.front_coordinate(), (coord_t{0ul, 0ul}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.front_coordinate()), coord_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).front_coordinate()), coord_t const &>));
    }
}

TYPED_TEST(alignment_result_test, alignment)
{
    using alignment_t = decltype(std::declval<TypeParam>().alignment);
    aligned_seq_type seq{'A'_dna4, 'T'_dna4, seqan3::gap{}, 'C'_dna4, seqan3::gap{}, seqan3::gap{}, 'A'_dna4};

    {
        seqan3::alignment_result<TypeParam> tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.alignment(), (alignment_t{seq, seq}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.alignment()), alignment_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).alignment()), alignment_t const &>));
    }

    {
        seqan3::alignment_result<TypeParam> const tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.alignment(), (alignment_t{seq, seq}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.alignment()), alignment_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).alignment()), alignment_t const &>));
    }

    if constexpr (seqan3::tuple_like<alignment_t>)
    {
        seqan3::alignment_result<TypeParam> tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(std::get<0>(tmp.alignment()) | seqan3::views::persist
                                               | seqan3::views::to_char
                                               | seqan3::views::to<std::string>,
                  std::string{"AT-C--A"});
        EXPECT_EQ(std::get<1>(tmp.alignment()) | seqan3::views::persist
                                               | seqan3::views::to_char
                                               | seqan3::views::to<std::string>,
                  std::string{"AT-C--A"});
    }
    else
    {
        seqan3::alignment_result<TypeParam> tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.alignment()[0] | seqan3::views::persist | seqan3::views::to_char | seqan3::views::to<std::string>,
                  "AT-C--A");
        EXPECT_EQ(tmp.alignment()[1] | seqan3::views::persist | seqan3::views::to_char | seqan3::views::to<std::string>,
                  "AT-C--A");
    }
}

TEST(alignment_result_test, reduced_type)
{
    {
        seqan3::detail::alignment_result_value_type tr{2u, 5};
        seqan3::alignment_result tmp(tr);
        EXPECT_EQ(tmp.id(), 2u);
        EXPECT_EQ(tmp.score(), 5);
    }

    {
        seqan3::detail::alignment_result_value_type tr{2, 5.0f, std::pair<int, int>{1, -1}};
        seqan3::alignment_result tmp(tr);
        EXPECT_EQ(tmp.id(), 2);
        EXPECT_FLOAT_EQ(tmp.score(), 5.0f);
        EXPECT_EQ((tmp.back_coordinate()), (std::pair<int, int>{1, -1}));
    }

    {
        seqan3::detail::alignment_result_value_type tr{2, 5.0f, std::pair<int, int>{1, -1}, std::pair<int, int>{10,
                                                                                                                -10}};
        seqan3::alignment_result tmp(tr);
        EXPECT_EQ(tmp.id(), 2);
        EXPECT_FLOAT_EQ(tmp.score(), 5.0f);
        EXPECT_EQ((tmp.back_coordinate()), (std::pair<int, int>{1, -1}));
        EXPECT_EQ((tmp.front_coordinate()), (std::pair<int, int>{10, -10}));
    }
}

TEST(alignment_result_test, type_deduction)
{
    {
        using coord_t = std::pair<int, int>;
        std::vector<seqan3::gapped<seqan3::rna5>> seq{'A'_rna5,
                                                      'U'_rna5,
                                                      seqan3::gap{},
                                                      'C'_rna5,
                                                      seqan3::gap{},
                                                      seqan3::gap{},
                                                      'A'_rna5};

        seqan3::detail::alignment_result_value_type tr{2, 5.0, coord_t{1, -1}, coord_t{10,-10}, seq};
        seqan3::alignment_result tmp(tr);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.id()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.score()), double>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.back_coordinate()), coord_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.front_coordinate()), coord_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.alignment()), std::vector<seqan3::gapped<seqan3::rna5>> const &>));

        EXPECT_EQ(tmp.id(), 2);
        EXPECT_DOUBLE_EQ(tmp.score(), 5.0);
        EXPECT_EQ((tmp.back_coordinate()), (coord_t{1, -1}));
        EXPECT_EQ((tmp.front_coordinate()), (coord_t{10, -10}));
        EXPECT_EQ(tmp.alignment(), seq);
    }

    {
        using coord_t = std::pair<int, int>;

        seqan3::detail::alignment_result_value_type tr{2, 5.0, coord_t{1, -1}};
        seqan3::alignment_result tmp(tr);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.id()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.score()), double>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.back_coordinate()), coord_t const &>));
    }

    {
        seqan3::detail::alignment_result_value_type tr{2, 5.0};
        seqan3::alignment_result tmp(tr);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.id()), int>));
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
    seqan3::detail::alignment_result_value_type result_value{2u, 5};
    seqan3::alignment_result result(result_value);

    using result_value_t = typename seqan3::detail::alignment_result_value_type_accessor<decltype(result)>::type;
    EXPECT_TRUE((std::same_as<result_value_t, decltype(result_value)>));
}
