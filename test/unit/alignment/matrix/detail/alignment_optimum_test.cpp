// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alignment/matrix/detail/alignment_optimum.hpp>
#include <seqan3/test/simd_utility.hpp>
#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/simd/simd.hpp>

template <typename scalar_t>
struct extract_scalar_type
{
    using type = scalar_t;
};

template <seqan3::simd::simd_concept simd_t>
struct extract_scalar_type<simd_t>
{
    using type = typename seqan3::simd::simd_traits<simd_t>::scalar_type;
};

template <typename test_t>
struct alignment_optimum_test : public ::testing::Test
{
    using scalar_t = typename extract_scalar_type<test_t>::type;

    template <typename lhs_t, typename rhs_t>
    void expect_eq(lhs_t lhs, rhs_t rhs)
    {
        if constexpr (seqan3::simd::simd_concept<lhs_t> && seqan3::simd::simd_concept<rhs_t>)
            SIMD_EQ(lhs, rhs);
        else if constexpr (seqan3::simd::simd_concept<lhs_t> && std::integral<rhs_t>)
            SIMD_EQ(lhs, seqan3::simd::fill<lhs_t>(rhs));
        else
            EXPECT_EQ(lhs, static_cast<lhs_t>(rhs));
    }
};

using score_types = ::testing::Types<int32_t, seqan3::simd::simd_type_t<int32_t>>;

TYPED_TEST_SUITE(alignment_optimum_test, score_types, );

TYPED_TEST(alignment_optimum_test, construction)
{
    using alignment_optimum_t = seqan3::detail::alignment_optimum<TypeParam>;

    EXPECT_TRUE(std::is_nothrow_default_constructible_v<TypeParam>);
    EXPECT_TRUE(std::is_nothrow_default_constructible_v<alignment_optimum_t>);
    EXPECT_TRUE(std::is_nothrow_copy_constructible_v<alignment_optimum_t>);
    EXPECT_TRUE(std::is_nothrow_copy_assignable_v<alignment_optimum_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<alignment_optimum_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<alignment_optimum_t>);
    EXPECT_TRUE(std::is_destructible_v<alignment_optimum_t>);
}

TYPED_TEST(alignment_optimum_test, type_deduction)
{
    seqan3::detail::alignment_optimum default_optimum{};
    EXPECT_TRUE((std::is_same_v<decltype(default_optimum), seqan3::detail::alignment_optimum<int32_t>>));

    seqan3::detail::alignment_optimum deduced_optimum{TypeParam{1}, TypeParam{2}, TypeParam{10}};
    EXPECT_TRUE((std::is_same_v<decltype(deduced_optimum), seqan3::detail::alignment_optimum<TypeParam>>));
}

TYPED_TEST(alignment_optimum_test, default_constructed)
{
    using scalar_t = typename TestFixture::scalar_t;
    seqan3::detail::alignment_optimum<TypeParam> default_optimum{};

    this->expect_eq(default_optimum.score, std::numeric_limits<scalar_t>::lowest());
    this->expect_eq(default_optimum.column_index, 0u);
    this->expect_eq(default_optimum.row_index, 0u);
}

TYPED_TEST(alignment_optimum_test, general_construction)
{
    seqan3::detail::alignment_optimum optimum{TypeParam{1}, TypeParam{2}, TypeParam{10}};

    this->expect_eq(optimum.score, TypeParam{10});
    this->expect_eq(optimum.column_index, TypeParam{1});
    this->expect_eq(optimum.row_index, TypeParam{2});
}

TYPED_TEST(alignment_optimum_test, update_if_new_optimal_score)
{
    using scalar_t = typename TestFixture::scalar_t;
    seqan3::detail::alignment_optimum<TypeParam> optimum{};

    this->expect_eq(optimum.score, std::numeric_limits<scalar_t>::lowest());
    this->expect_eq(optimum.column_index, 0u);
    this->expect_eq(optimum.row_index, 0u);

    // Bigger score.
    optimum.update_if_new_optimal_score(TypeParam{10},
                                        seqan3::detail::column_index_type{1},
                                        seqan3::detail::row_index_type{2});

    this->expect_eq(optimum.score, TypeParam{10});
    this->expect_eq(optimum.column_index, 1u);
    this->expect_eq(optimum.row_index, 2u);

    // Same score.
    optimum.update_if_new_optimal_score(TypeParam{10},
                                        seqan3::detail::column_index_type{4},
                                        seqan3::detail::row_index_type{5});

    this->expect_eq(optimum.score, TypeParam{10});
    this->expect_eq(optimum.column_index, 1u);
    this->expect_eq(optimum.row_index, 2u);

    // Lower score.
    optimum.update_if_new_optimal_score(TypeParam{7},
                                        seqan3::detail::column_index_type{4},
                                        seqan3::detail::row_index_type{5});

    this->expect_eq(optimum.score, TypeParam{10});
    this->expect_eq(optimum.column_index, 1u);
    this->expect_eq(optimum.row_index, 2u);

    // Mixed score differences
    if constexpr (seqan3::simd::simd_concept<TypeParam>)
    { // The following will only work if the simd type has more than one element.
        if constexpr (seqan3::simd::simd_traits<TypeParam>::length > 1)
        {
            TypeParam score_vector{5};
            TypeParam cmp_col_index = optimum.column_index;
            TypeParam cmp_row_index = optimum.row_index;

            score_vector[1] = 11;
            cmp_col_index[1] = 3;
            cmp_row_index[1] = 7;

            optimum.update_if_new_optimal_score(score_vector,
                                                seqan3::detail::column_index_type{3},
                                                seqan3::detail::row_index_type{7});

            TypeParam cmp_score_vector{10};
            cmp_score_vector[1] = 11;
            this->expect_eq(optimum.score, cmp_score_vector);
            this->expect_eq(optimum.column_index, cmp_col_index);
            this->expect_eq(optimum.row_index, cmp_row_index);
        }
    }
}
