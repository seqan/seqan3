// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/detail/simd_matrix_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/test/pretty_printing.hpp>
#include <seqan3/test/simd_utility.hpp>

template <typename simd_t>
struct simd_matrix_scoring_scheme_test : public ::testing::Test
{
    using scalar_t = typename seqan3::simd_traits<simd_t>::scalar_type;

    scalar_t padded_value1 = std::numeric_limits<scalar_t>::lowest();      // sets the most significant bit.
    scalar_t padded_value2 = std::numeric_limits<scalar_t>::lowest() >> 1; // sets the bit before most significant bit.
};

using simd_test_types = ::testing::Types<//seqan3::simd::simd_type_t<int8_t>,
                                         //seqan3::simd::simd_type_t<int16_t>,
                                         seqan3::simd::simd_type_t<int32_t>>;

TYPED_TEST_SUITE(simd_matrix_scoring_scheme_test, simd_test_types, );

TYPED_TEST(simd_matrix_scoring_scheme_test, basic_construction)
{
    using scheme_t = seqan3::detail::simd_matrix_scoring_scheme<TypeParam,
                                                                seqan3::aa27,
                                                                seqan3::detail::global_alignment_type,
                                                                seqan3::aminoacid_scoring_scheme<>>;

    EXPECT_TRUE(std::is_nothrow_default_constructible_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_copy_constructible_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_copy_assignable_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<scheme_t>);
    EXPECT_TRUE((std::is_constructible_v<scheme_t, seqan3::aminoacid_scoring_scheme<>>));
    EXPECT_TRUE(std::semiregular<scheme_t>);
}

TYPED_TEST(simd_matrix_scoring_scheme_test, construct_from_scoring_scheme_nothrow)
{
    using scheme_t = seqan3::detail::simd_matrix_scoring_scheme<TypeParam,
                                                                seqan3::aa27,
                                                                seqan3::detail::global_alignment_type,
                                                                seqan3::aminoacid_scoring_scheme<>>;

    scheme_t simd_scheme{seqan3::aminoacid_scoring_scheme{seqan3::aminoacid_similarity_matrix::BLOSUM30}};

    TypeParam simd_value1 = seqan3::simd::fill<TypeParam>(2);
    TypeParam simd_value2 = seqan3::simd::fill<TypeParam>(2);
    SIMD_EQ(simd_scheme.score(simd_value1, simd_value2), seqan3::simd::fill<TypeParam>(17));

    simd_value2 = seqan3::simd::fill<TypeParam>(1);
    SIMD_EQ(simd_scheme.score(simd_value1, simd_value2), seqan3::simd::fill<TypeParam>(-2));
}

TYPED_TEST(simd_matrix_scoring_scheme_test, construct_from_scoring_scheme_throw_on_overflow)
{
    using scalar_t = typename seqan3::simd_traits<TypeParam>::scalar_type;
    using scheme_t = seqan3::detail::simd_matrix_scoring_scheme<TypeParam,
                                                                seqan3::aa27,
                                                                seqan3::detail::global_alignment_type,
                                                                seqan3::aminoacid_scoring_scheme<int64_t>>;

    int64_t too_big = static_cast<int64_t>(std::numeric_limits<scalar_t>::max()) + 1;
    int64_t too_small = static_cast<int64_t>(std::numeric_limits<scalar_t>::lowest()) - 1;

    typename seqan3::aminoacid_scoring_scheme<int64_t>::matrix_type matrix{};

    EXPECT_NO_THROW(scheme_t{seqan3::aminoacid_scoring_scheme<int64_t>{matrix}});

    matrix[0][0] = too_big;
    EXPECT_THROW((scheme_t{seqan3::aminoacid_scoring_scheme<int64_t>{matrix}}),
                 std::invalid_argument);
    matrix[0][0] = too_small;
    EXPECT_THROW((scheme_t{seqan3::aminoacid_scoring_scheme<int64_t>{matrix}}),
                 std::invalid_argument);

    matrix[0][0] = 0;
    matrix[seqan3::alphabet_size<seqan3::aa27> - 1][seqan3::alphabet_size<seqan3::aa27> - 1] = too_big;
    EXPECT_THROW((scheme_t{seqan3::aminoacid_scoring_scheme<int64_t>{matrix}}),
                 std::invalid_argument);
    matrix[seqan3::alphabet_size<seqan3::aa27> - 1][seqan3::alphabet_size<seqan3::aa27> - 1] = too_small;
    EXPECT_THROW((scheme_t{seqan3::aminoacid_scoring_scheme<int64_t>{matrix}}),
                 std::invalid_argument);
}

TYPED_TEST(simd_matrix_scoring_scheme_test, score_global)
{
    using scheme_t = seqan3::detail::simd_matrix_scoring_scheme<TypeParam,
                                                                seqan3::aa27,
                                                                seqan3::detail::global_alignment_type,
                                                                seqan3::aminoacid_scoring_scheme<>>;

    scheme_t scheme{seqan3::aminoacid_scoring_scheme{seqan3::aminoacid_similarity_matrix::BLOSUM30}};

    TypeParam simd_value1 = seqan3::simd::fill<TypeParam>(2);
    TypeParam simd_value2 = seqan3::simd::fill<TypeParam>(2);

    // all match
    SIMD_EQ(scheme.score(simd_value1, simd_value2), seqan3::simd::fill<TypeParam>(17));

    // all mismatch
    simd_value2 = seqan3::simd::fill<TypeParam>(3);
    SIMD_EQ(scheme.score(simd_value1, simd_value2), seqan3::simd::fill<TypeParam>(-3));

    // first matches, remaining mismatches.
    simd_value2[0] = 2;
    TypeParam result = seqan3::simd::fill<TypeParam>(-3);
    result[0] = 17;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // first mismatches, remaining matches.
    simd_value1 = simd_value2;
    simd_value1[0] = 3;
    result = seqan3::simd::fill<TypeParam>(9);
    result[0] = -3;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);
}

TYPED_TEST(simd_matrix_scoring_scheme_test, score_global_with_padding)
{
    using scheme_t = seqan3::detail::simd_matrix_scoring_scheme<TypeParam,
                                                                seqan3::aa27,
                                                                seqan3::detail::global_alignment_type,
                                                                seqan3::aminoacid_scoring_scheme<>>;

    scheme_t scheme{seqan3::aminoacid_scoring_scheme{seqan3::aminoacid_similarity_matrix::BLOSUM30}};

    TypeParam simd_value1 = seqan3::simd::fill<TypeParam>(2);
    TypeParam simd_value2 = seqan3::simd::fill<TypeParam>(3);
    TypeParam result = seqan3::simd::fill<TypeParam>(-3);

    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // First value is regular symbol; second value is padded symbol => score of 1.
    simd_value2[0] = this->padded_value1;
    result[0] = 1;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // First value is padded symbol; second value is padded symbol => score of 1.
    simd_value1[0] = this->padded_value1;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // First value is padded symbol; second value is regular symbol => score of 1.
    simd_value2[0] = 2;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);
}

TYPED_TEST(simd_matrix_scoring_scheme_test, score_local)
{
    // In local alignment we always want to mismatch.
    using scheme_t = seqan3::detail::simd_matrix_scoring_scheme<TypeParam,
                                                                seqan3::aa27,
                                                                seqan3::detail::local_alignment_type,
                                                                seqan3::aminoacid_scoring_scheme<>>;

    scheme_t scheme{seqan3::aminoacid_scoring_scheme{seqan3::aminoacid_similarity_matrix::BLOSUM30}};

    TypeParam simd_value1 = seqan3::simd::fill<TypeParam>(2);
    TypeParam simd_value2 = seqan3::simd::fill<TypeParam>(2);

    // all match
    SIMD_EQ(scheme.score(simd_value1, simd_value2), seqan3::simd::fill<TypeParam>(17));

    // all mismatch
    simd_value2 = seqan3::simd::fill<TypeParam>(3);
    SIMD_EQ(scheme.score(simd_value1, simd_value2), seqan3::simd::fill<TypeParam>(-3));

    // first matches, remaining mismatches.
    simd_value2[0] = 2;
    TypeParam result = seqan3::simd::fill<TypeParam>(-3);
    result[0] = 17;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // first mismatches, remaining matches.
    simd_value1 = simd_value2;
    simd_value1[0] = 3;
    result = seqan3::simd::fill<TypeParam>(9);
    result[0] = -3;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);
}

TYPED_TEST(simd_matrix_scoring_scheme_test, score_local_with_padding)
{
    // In local alignment we always want to mismatch.
    using scheme_t = seqan3::detail::simd_matrix_scoring_scheme<TypeParam,
                                                                seqan3::aa27,
                                                                seqan3::detail::local_alignment_type,
                                                                seqan3::aminoacid_scoring_scheme<>>;

    scheme_t scheme{seqan3::aminoacid_scoring_scheme{seqan3::aminoacid_similarity_matrix::BLOSUM30}};

    TypeParam simd_value1 = seqan3::simd::fill<TypeParam>(2);
    TypeParam simd_value2 = seqan3::simd::fill<TypeParam>(2);
    TypeParam result = seqan3::simd::fill<TypeParam>(17);

    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // First value is regular symbol; second value is padded symbol => score of -1.
    simd_value2[0] = this->padded_value2;
    result[0] = -1;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // First value is padded symbol; second value is padded symbol => score of -1.
    simd_value1[0] = this->padded_value1;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // First value is padded symbol; second value is regular symbol => score of -1.
    simd_value2[0] = 3;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);
}
