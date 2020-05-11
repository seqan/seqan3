// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/detail/simd_match_mismatch_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/test/pretty_printing.hpp>
#include <seqan3/test/simd_utility.hpp>

template <typename simd_t>
struct simd_match_mismatch_scoring_scheme_test : public ::testing::Test
{
    using scalar_t = typename seqan3::simd_traits<simd_t>::scalar_type;

    scalar_t padded_value1 = std::numeric_limits<scalar_t>::lowest();      // sets the most significant bit.
    scalar_t padded_value2 = std::numeric_limits<scalar_t>::lowest() >> 1; // sets the bit before most significant bit.
};

using simd_test_types = ::testing::Types<seqan3::simd::simd_type_t<int8_t>,
                                         seqan3::simd::simd_type_t<int16_t>,
                                         seqan3::simd::simd_type_t<int32_t>>;

TYPED_TEST_SUITE(simd_match_mismatch_scoring_scheme_test, simd_test_types, );

TYPED_TEST(simd_match_mismatch_scoring_scheme_test, basic_construction)
{
    using scheme_t = seqan3::detail::simd_match_mismatch_scoring_scheme<TypeParam,
                                                                        seqan3::dna4,
                                                                        seqan3::detail::global_alignment_type>;

    EXPECT_TRUE(std::is_nothrow_default_constructible_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_copy_constructible_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_copy_assignable_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<scheme_t>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<scheme_t>);
    EXPECT_TRUE((std::is_constructible_v<scheme_t, seqan3::nucleotide_scoring_scheme<int16_t>>));
    EXPECT_TRUE(std::semiregular<scheme_t>);
}

TYPED_TEST(simd_match_mismatch_scoring_scheme_test, construct_from_scoring_scheme_nothrow)
{
    using scheme_t = seqan3::detail::simd_match_mismatch_scoring_scheme<TypeParam,
                                                                        seqan3::dna4,
                                                                        seqan3::detail::global_alignment_type>;

    scheme_t simd_scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4}, seqan3::mismatch_score{-5}}};

    TypeParam simd_value1 = seqan3::simd::fill<TypeParam>(2);
    TypeParam simd_value2 = seqan3::simd::fill<TypeParam>(2);
    SIMD_EQ(simd_scheme.score(simd_value1, simd_value2), seqan3::simd::fill<TypeParam>(4));

    simd_value2 = seqan3::simd::fill<TypeParam>(1);
    SIMD_EQ(simd_scheme.score(simd_value1, simd_value2), seqan3::simd::fill<TypeParam>(-5));
}

TYPED_TEST(simd_match_mismatch_scoring_scheme_test, construct_from_scoring_scheme_throw_on_overflow)
{
    using scalar_t = typename seqan3::simd_traits<TypeParam>::scalar_type;
    using scheme_t = seqan3::detail::simd_match_mismatch_scoring_scheme<TypeParam,
                                                                        seqan3::dna4,
                                                                        seqan3::detail::global_alignment_type>;

    int64_t too_big = static_cast<int64_t>(std::numeric_limits<scalar_t>::max()) + 1;
    int64_t too_small = static_cast<int64_t>(std::numeric_limits<scalar_t>::lowest()) - 1;
    EXPECT_THROW((scheme_t{seqan3::nucleotide_scoring_scheme<int64_t>{seqan3::match_score{too_big},
                                                                      seqan3::mismatch_score<int64_t>{-5}}}),
                 std::invalid_argument);
    EXPECT_THROW((scheme_t{seqan3::nucleotide_scoring_scheme<int64_t>{seqan3::match_score<int64_t>{4},
                                                                      seqan3::mismatch_score{too_small}}}),
                 std::invalid_argument);
}

TYPED_TEST(simd_match_mismatch_scoring_scheme_test, score_global)
{
    using scheme_t = seqan3::detail::simd_match_mismatch_scoring_scheme<TypeParam,
                                                                        seqan3::dna4,
                                                                        seqan3::detail::global_alignment_type>;

    scheme_t scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4}, seqan3::mismatch_score{-5}}};

    TypeParam simd_value1 = seqan3::simd::fill<TypeParam>(2);
    TypeParam simd_value2 = seqan3::simd::fill<TypeParam>(2);

    // all match
    SIMD_EQ(scheme.score(simd_value1, simd_value2), seqan3::simd::fill<TypeParam>(4));

    // all mismatch
    simd_value2 = seqan3::simd::fill<TypeParam>(3);
    SIMD_EQ(scheme.score(simd_value1, simd_value2), seqan3::simd::fill<TypeParam>(-5));

    // first matches, remaining mismatches.
    simd_value2[0] = 2;
    TypeParam result = seqan3::simd::fill<TypeParam>(-5);
    result[0] = 4;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // first mismatches, remaining matches.
    simd_value1 = simd_value2;
    simd_value1[0] = 1;
    result = seqan3::simd::fill<TypeParam>(4);
    result[0] = -5;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);
}

TYPED_TEST(simd_match_mismatch_scoring_scheme_test, score_global_with_padding)
{
    using scheme_t = seqan3::detail::simd_match_mismatch_scoring_scheme<TypeParam,
                                                                        seqan3::dna4,
                                                                        seqan3::detail::global_alignment_type>;

    scheme_t scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4}, seqan3::mismatch_score{-5}}};

    TypeParam simd_value1 = seqan3::simd::fill<TypeParam>(2);
    TypeParam simd_value2 = seqan3::simd::fill<TypeParam>(1);
    TypeParam result = seqan3::simd::fill<TypeParam>(-5);

    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // First value is regular symbol; second value is padded symbol => match.
    simd_value2[0] = this->padded_value1;
    result[0] = 4;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // First value is padded symbol; second value is padded symbol => match.
    simd_value1[0] = this->padded_value1;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // First value is padded symbol; second value is regular symbol => match.
    simd_value2[0] = 3;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);
}

TYPED_TEST(simd_match_mismatch_scoring_scheme_test, score_local)
{
    // In local alignment we always want to mismatch.
    using scheme_t = seqan3::detail::simd_match_mismatch_scoring_scheme<TypeParam,
                                                                        seqan3::dna4,
                                                                        seqan3::detail::local_alignment_type>;

    scheme_t scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4}, seqan3::mismatch_score{-5}}};

    TypeParam simd_value1 = seqan3::simd::fill<TypeParam>(2);
    TypeParam simd_value2 = seqan3::simd::fill<TypeParam>(2);

    // all match
    SIMD_EQ(scheme.score(simd_value1, simd_value2), seqan3::simd::fill<TypeParam>(4));

    // all mismatch
    simd_value2 = seqan3::simd::fill<TypeParam>(3);
    SIMD_EQ(scheme.score(simd_value1, simd_value2), seqan3::simd::fill<TypeParam>(-5));

    // first matches, remaining mismatches.
    simd_value2[0] = 2;
    TypeParam result = seqan3::simd::fill<TypeParam>(-5);
    result[0] = 4;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // first mismatches, remaining matches.
    simd_value1 = simd_value2;
    simd_value1[0] = 1;
    result = seqan3::simd::fill<TypeParam>(4);
    result[0] = -5;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);
}

TYPED_TEST(simd_match_mismatch_scoring_scheme_test, score_local_with_padding)
{
    // In local alignment we always want to mismatch.
    using scheme_t = seqan3::detail::simd_match_mismatch_scoring_scheme<TypeParam,
                                                                        seqan3::dna4,
                                                                        seqan3::detail::local_alignment_type>;

    scheme_t scheme{seqan3::nucleotide_scoring_scheme{seqan3::match_score{4}, seqan3::mismatch_score{-5}}};

    TypeParam simd_value1 = seqan3::simd::fill<TypeParam>(2);
    TypeParam simd_value2 = seqan3::simd::fill<TypeParam>(2);
    TypeParam result = seqan3::simd::fill<TypeParam>(4);

    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // First value is regular symbol; second value is padded symbol => mismatch.
    simd_value2[0] = this->padded_value2;
    result[0] = -5;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // First value is padded symbol; second value is padded symbol => mismatch.
    simd_value1[0] = this->padded_value1;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);

    // First value is padded symbol; second value is regular symbol => mismatch.
    simd_value2[0] = 3;
    SIMD_EQ(scheme.score(simd_value1, simd_value2), result);
}
