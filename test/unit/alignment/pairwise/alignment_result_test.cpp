// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
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
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/view/to_char.hpp>

using namespace seqan3;

using aligned_seq_type = std::vector<gapped<dna4>>;

template <typename T>
struct alignment_result_test : public ::testing::Test
{};

using alignment_result_test_types = ::testing::Types
      <detail::alignment_result_value_type<uint32_t, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>,
                                       std::pair<aligned_seq_type, aligned_seq_type>>,
       detail::alignment_result_value_type<uint32_t, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>,
                                       std::tuple<aligned_seq_type, aligned_seq_type>>,
       detail::alignment_result_value_type<uint32_t, int32_t, std::pair<size_t, size_t>, std::pair<size_t, size_t>,
                                       std::vector<aligned_seq_type>>,
       detail::alignment_result_value_type<uint32_t, float,   std::pair<size_t, size_t>, std::pair<size_t, size_t>,
                                       std::pair<aligned_seq_type, aligned_seq_type>>,
       detail::alignment_result_value_type<uint32_t, float,   std::pair<size_t, size_t>, std::pair<size_t, size_t>,
                                       std::tuple<aligned_seq_type, aligned_seq_type>>,
       detail::alignment_result_value_type<uint32_t, float,   std::pair<size_t, size_t>, std::pair<size_t, size_t>,
                                       std::vector<aligned_seq_type>>>;

TYPED_TEST_CASE(alignment_result_test, alignment_result_test_types);

TYPED_TEST(alignment_result_test, type_specialisation)
{
    EXPECT_TRUE((detail::is_type_specialisation_of_v<TypeParam, detail::alignment_result_value_type>));
}

TYPED_TEST(alignment_result_test, constructor)
{
    EXPECT_TRUE((std::is_default_constructible_v<alignment_result<TypeParam>>));
    EXPECT_TRUE((std::is_copy_constructible_v<alignment_result<TypeParam>>));
    EXPECT_TRUE((std::is_move_constructible_v<alignment_result<TypeParam>>));
    EXPECT_TRUE((std::is_copy_assignable_v<alignment_result<TypeParam>>));
    EXPECT_TRUE((std::is_move_assignable_v<alignment_result<TypeParam>>));
    EXPECT_TRUE((std::is_destructible_v<alignment_result<TypeParam>>));
}

TYPED_TEST(alignment_result_test, get_id)
{
    using id_t = decltype(std::declval<TypeParam>().id);
    aligned_seq_type seq{'A'_dna4, 'T'_dna4, gap{}, 'C'_dna4, gap{}, gap{}, 'A'_dna4};

    {
        alignment_result<TypeParam> tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_TRUE((std::is_same_v<decltype(tmp.id()), id_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).id()), id_t>));
        EXPECT_EQ(tmp.id(), 1u);
    }

    {
        alignment_result<TypeParam> const tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.id(), 1u);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.id()), id_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).id()), id_t>));
    }
}

TYPED_TEST(alignment_result_test, get_score)
{
    using score_t = decltype(std::declval<TypeParam>().score);
    aligned_seq_type seq{'A'_dna4, 'T'_dna4, gap{}, 'C'_dna4, gap{}, gap{}, 'A'_dna4};

    {
        alignment_result<TypeParam> tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.score(), 0);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.score()), score_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).score()), score_t>));
    }

    {
        alignment_result<TypeParam> const tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.score(), 0);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.score()), score_t>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).score()), score_t>));
    }

}

TYPED_TEST(alignment_result_test, back_coordinate)
{
    using coord_t = decltype(std::declval<TypeParam>().back_coordinate);
    aligned_seq_type seq{'A'_dna4, 'T'_dna4, gap{}, 'C'_dna4, gap{}, gap{}, 'A'_dna4};

    {
        alignment_result<TypeParam> tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.back_coordinate(), (coord_t{10ul, 10ul}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.back_coordinate()), coord_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).back_coordinate()), coord_t const &>));
    }

    {
        alignment_result<TypeParam> const tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.back_coordinate(), (coord_t{10ul, 10ul}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.back_coordinate()), coord_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).back_coordinate()), coord_t const &>));
    }
}

TYPED_TEST(alignment_result_test, front_coordinate)
{
    using coord_t = decltype(std::declval<TypeParam>().front_coordinate);
    aligned_seq_type seq{'A'_dna4, 'T'_dna4, gap{}, 'C'_dna4, gap{}, gap{}, 'A'_dna4};

    {
        alignment_result<TypeParam> tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.front_coordinate(), (coord_t{0ul, 0ul}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.front_coordinate()), coord_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).front_coordinate()), coord_t const &>));
    }

    {
        alignment_result<TypeParam> const tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.front_coordinate(), (coord_t{0ul, 0ul}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.front_coordinate()), coord_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).front_coordinate()), coord_t const &>));
    }
}

TYPED_TEST(alignment_result_test, alignment)
{
    using alignment_t = decltype(std::declval<TypeParam>().alignment);
    aligned_seq_type seq{'A'_dna4, 'T'_dna4, gap{}, 'C'_dna4, gap{}, gap{}, 'A'_dna4};

    {
        alignment_result<TypeParam> tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.alignment(), (alignment_t{seq, seq}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.alignment()), alignment_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).alignment()), alignment_t const &>));
    }

    {
        alignment_result<TypeParam> const tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(tmp.alignment(), (alignment_t{seq, seq}));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.alignment()), alignment_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(std::move(tmp).alignment()), alignment_t const &>));
    }

    if constexpr (TupleLike<alignment_t>)
    {
        alignment_result<TypeParam> tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(std::string{std::get<0>(tmp.alignment()) | view::persist | view::to_char},
                  std::string{"AT-C--A"});
        EXPECT_EQ(std::string{std::get<1>(tmp.alignment()) | view::persist | view::to_char},
                  std::string{"AT-C--A"});
    }
    else
    {
        alignment_result<TypeParam> tmp{TypeParam{1u, 0, {10ul, 10ul}, {0ul, 0ul}, {seq, seq}}};
        EXPECT_EQ(std::string{(tmp.alignment()[0]) | view::persist | view::to_char}, std::string{"AT-C--A"});
        EXPECT_EQ(std::string{(tmp.alignment()[1]) | view::persist | view::to_char}, std::string{"AT-C--A"});
    }
}

TEST(alignment_result_test, reduced_type)
{
    {
        detail::alignment_result_value_type tr{2u, 5};
        alignment_result tmp(tr);
        EXPECT_EQ(tmp.id(), 2u);
        EXPECT_EQ(tmp.score(), 5);
    }

    {
        detail::alignment_result_value_type tr{2, 5.0f, std::pair<int, int>{1, -1}};
        alignment_result tmp(tr);
        EXPECT_EQ(tmp.id(), 2);
        EXPECT_FLOAT_EQ(tmp.score(), 5.0f);
        EXPECT_EQ((tmp.back_coordinate()), (std::pair<int, int>{1, -1}));
    }

    {
        detail::alignment_result_value_type tr{2, 5.0f, std::pair<int, int>{1, -1}, std::pair<int, int>{10, -10}};
        alignment_result tmp(tr);
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
        std::vector<gapped<rna5>> seq{'A'_rna5, 'U'_rna5, gap{}, 'C'_rna5, gap{}, gap{}, 'A'_rna5};

        detail::alignment_result_value_type tr{2, 5.0, coord_t{1, -1}, coord_t{10,-10}, seq};
        alignment_result tmp(tr);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.id()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.score()), double>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.back_coordinate()), coord_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.front_coordinate()), coord_t const &>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.alignment()), std::vector<gapped<rna5>> const &>));

        EXPECT_EQ(tmp.id(), 2);
        EXPECT_DOUBLE_EQ(tmp.score(), 5.0);
        EXPECT_EQ((tmp.back_coordinate()), (coord_t{1, -1}));
        EXPECT_EQ((tmp.front_coordinate()), (coord_t{10, -10}));
        EXPECT_EQ(tmp.alignment(), seq);
    }

    {
        using coord_t = std::pair<int, int>;

        detail::alignment_result_value_type tr{2, 5.0, coord_t{1, -1}};
        alignment_result tmp(tr);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.id()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.score()), double>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.back_coordinate()), coord_t const &>));
    }

    {
        detail::alignment_result_value_type tr{2, 5.0};
        alignment_result tmp(tr);
        EXPECT_TRUE((std::is_same_v<decltype(tmp.id()), int>));
        EXPECT_TRUE((std::is_same_v<decltype(tmp.score()), double>));
    }
}

TEST(alignment_result_test, empty_type)
{
    detail::alignment_result_value_type tr{};
    alignment_result tmp(tr);
    // members must not be accessed
}
