// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/scoring/simd_scoring_scheme_base.hpp>
#include <seqan3/test/pretty_printing.hpp>
#include <seqan3/test/simd_utility.hpp>

using namespace seqan3;

TEST(simd_scoring_scheme_simple_test, construction)
{
    using simd_t = typename simd::simd_type<int32_t>::type;
    using scheme_t = simd_scoring_scheme_simple<simd_t, detail::global_alignment_type>;

    EXPECT_TRUE(std::is_default_constructible_v<scheme_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<scheme_t>);
    EXPECT_TRUE(std::is_move_constructible_v<scheme_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<scheme_t>);
    EXPECT_TRUE(std::is_move_assignable_v<scheme_t>);
    EXPECT_TRUE(std::is_destructible_v<scheme_t>);
    EXPECT_TRUE((std::is_constructible_v<scheme_t, match_score<int32_t>, mismatch_score<int32_t>>));
    EXPECT_TRUE((std::is_constructible_v<scheme_t, match_score<int8_t>, mismatch_score<int8_t>>));
    EXPECT_TRUE((std::is_constructible_v<scheme_t, match_score<float>, mismatch_score<float>>));
    EXPECT_TRUE(std::Semiregular<scheme_t>);
}

TEST(simd_scoring_scheme_simple_test, set_simple_scheme)
{
    using simd_t = typename simd::simd_type<int32_t>::type;
    using scheme_t = simd_scoring_scheme_simple<simd_t, detail::global_alignment_type>;

    scheme_t scheme{};
    scheme.set_simple_scheme(match_score{4}, mismatch_score{-5});

    simd_t s1 = simd::fill<simd_t>(2);
    simd_t s2 = simd::fill<simd_t>(2);

    SIMD_EQ(scheme.score(s1, s2), simd::fill<simd_t>(4));

    s2 = simd::fill<simd_t>(3);
    SIMD_EQ(scheme.score(s1, s2), simd::fill<simd_t>(-5));
}

TEST(simd_scoring_scheme_simple_test, score_global)
{
    using simd_t = typename simd::simd_type<int32_t>::type;
    using scheme_t = simd_scoring_scheme_simple<simd_t, detail::global_alignment_type>;

    scheme_t scheme{match_score{4}, mismatch_score{-5}};

    simd_t s1 = simd::fill<simd_t>(2);
    simd_t s2 = simd::fill<simd_t>(2);

    // all match
    SIMD_EQ(scheme.score(s1, s2), simd::fill<simd_t>(4));

    // all mismatch
    s2 = simd::fill<simd_t>(3);
    SIMD_EQ(scheme.score(s1, s2), simd::fill<simd_t>(-5));

    // first one is padded so it must be a match.
    s2[0] = 1 << (sizeof(typename simd_traits<simd_t>::scalar_type) * 8 - 1);
    simd_t res = simd::fill<simd_t>(-5);
    res[0] = 4;
    SIMD_EQ(scheme.score(s1, s2), res);

    // first one in other sequence is padded so it must be still a match.
    s1[0] = 1 << (sizeof(typename simd_traits<simd_t>::scalar_type) * 8 - 1);
    SIMD_EQ(scheme.score(s1, s2), res);

    // Only first one in other sequence is padded so it must be still a match.
    s2[0] = 3;
    SIMD_EQ(scheme.score(s1, s2), res);
}

TEST(simd_scoring_scheme_simple_test, score_local)
{
    //TODO Add as soon local alignment is in the library.
}
