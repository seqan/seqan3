// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/test/cereal.hpp>

TEST(gap_scheme, constructors_and_type_deduction_guides)
{

    {
        seqan3::gap_scheme scheme;
        EXPECT_TRUE((std::is_same_v<decltype(scheme), seqan3::gap_scheme<int8_t>>));
    }

    {
        seqan3::gap_scheme scheme{};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), seqan3::gap_scheme<int8_t>>));
    }

    {
        seqan3::gap_scheme scheme{seqan3::gap_score{-2}, seqan3::gap_open_score{-4}};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), seqan3::gap_scheme<int8_t>>));
    }

    {
        seqan3::gap_scheme scheme{seqan3::gap_score{-2}};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), seqan3::gap_scheme<int8_t>>));
    }

    {
        seqan3::gap_scheme scheme{seqan3::gap_score{-2.}, seqan3::gap_open_score{-4.}};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), seqan3::gap_scheme<float>>));
    }

    {
        seqan3::gap_scheme scheme{seqan3::gap_score{-2.}};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), seqan3::gap_scheme<float>>));
    }
}

TEST(gap_scheme, member_types)
{
    seqan3::gap_scheme scheme{};

    using score_t = typename decltype(scheme)::score_type;
    EXPECT_TRUE((std::is_same_v<score_t, int8_t>));
}

TEST(gap_scheme, get_gap_score)
{
    seqan3::gap_scheme scheme{};
    EXPECT_EQ(scheme.get_gap_score(), -1);
    EXPECT_TRUE((std::is_same_v<typename decltype(scheme)::score_type &, decltype(scheme.get_gap_score())>));
}

TEST(gap_scheme, set_score_gap)
{
    seqan3::gap_scheme scheme{};
    EXPECT_EQ(scheme.get_gap_score(), -1);
    scheme.get_gap_score() = -2;
    EXPECT_EQ(scheme.get_gap_score(), -2);
}

TEST(gap_scheme, get_gap_open_score)
{
    seqan3::gap_scheme scheme{};
    EXPECT_EQ(scheme.get_gap_open_score(), 0);
    EXPECT_TRUE((std::is_same_v<typename decltype(scheme)::score_type &, decltype(scheme.get_gap_open_score())>));
}

TEST(gap_scheme, set_score_gap_open)
{
    seqan3::gap_scheme scheme{};
    EXPECT_EQ(scheme.get_gap_open_score(), 0);
    scheme.get_gap_open_score() = -2;
    EXPECT_EQ(scheme.get_gap_open_score(), -2);
}

TEST(gap_scheme, set_linear)
{
    seqan3::gap_scheme scheme{seqan3::gap_score{-2}};
    EXPECT_EQ(scheme.get_gap_score(), -2);
    EXPECT_EQ(scheme.get_gap_open_score(), 0);

    scheme.set_linear(seqan3::gap_score{-3});
    EXPECT_EQ(scheme.get_gap_score(), -3);
    EXPECT_EQ(scheme.get_gap_open_score(), 0);
}

TEST(gap_scheme, set_affine)
{
    seqan3::gap_scheme scheme{seqan3::gap_score{-2}, seqan3::gap_open_score{-4}};
    EXPECT_EQ(scheme.get_gap_score(), -2);
    EXPECT_EQ(scheme.get_gap_open_score(), -4);

    scheme.set_affine(seqan3::gap_score{-3}, seqan3::gap_open_score{-6});
    EXPECT_EQ(scheme.get_gap_score(), -3);
    EXPECT_EQ(scheme.get_gap_open_score(), -6);
}

TEST(gap_scheme, score)
{
    seqan3::gap_scheme scheme{seqan3::gap_score{-2}};
    EXPECT_EQ(scheme.score(0), 0);
    EXPECT_EQ(scheme.score(2), -4);
    EXPECT_EQ(scheme.score(5), -10);

    scheme.set_affine(seqan3::gap_score{-3}, seqan3::gap_open_score{-6});
    EXPECT_EQ(scheme.score(0), 0);
    EXPECT_EQ(scheme.score(2), -12);
    EXPECT_EQ(scheme.score(5), -21);
}

TEST(gap_scheme, serialisation)
{
    seqan3::gap_scheme scheme1;

    scheme1.set_linear(seqan3::gap_score{-3});
    seqan3::test::do_serialisation(scheme1);

    scheme1.set_affine(seqan3::gap_score{-3}, seqan3::gap_open_score{-6});
    seqan3::test::do_serialisation(scheme1);
}
