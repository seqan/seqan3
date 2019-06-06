// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/test/cereal.hpp>

using namespace seqan3;

TEST(gap_scheme, constructors_and_type_deduction_guides)
{

    {
        gap_scheme scheme;
        EXPECT_TRUE((std::is_same_v<decltype(scheme), gap_scheme<int8_t>>));
    }

    {
        gap_scheme scheme{};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), gap_scheme<int8_t>>));
    }

    {
        gap_scheme scheme{gap_score{-2}, gap_open_score{-4}};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), gap_scheme<int8_t>>));
    }

    {
        gap_scheme scheme{gap_score{-2}};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), gap_scheme<int8_t>>));
    }

    {
        gap_scheme scheme{gap_score{-2.}, gap_open_score{-4.}};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), gap_scheme<float>>));
    }

    {
        gap_scheme scheme{gap_score{-2.}};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), gap_scheme<float>>));
    }
}

TEST(gap_scheme, member_types)
{
    gap_scheme scheme{};

    using score_t = typename decltype(scheme)::score_type;
    EXPECT_TRUE((std::is_same_v<score_t, int8_t>));
}

TEST(gap_scheme, get_gap_score)
{
    gap_scheme scheme{};
    EXPECT_EQ(scheme.get_gap_score(), -1);
    EXPECT_TRUE((std::is_same_v<typename decltype(scheme)::score_type &, decltype(scheme.get_gap_score())>));
}

TEST(gap_scheme, set_score_gap)
{
    gap_scheme scheme{};
    EXPECT_EQ(scheme.get_gap_score(), -1);
    scheme.get_gap_score() = -2;
    EXPECT_EQ(scheme.get_gap_score(), -2);
}

TEST(gap_scheme, get_gap_open_score)
{
    gap_scheme scheme{};
    EXPECT_EQ(scheme.get_gap_open_score(), 0);
    EXPECT_TRUE((std::is_same_v<typename decltype(scheme)::score_type &, decltype(scheme.get_gap_open_score())>));
}

TEST(gap_scheme, set_score_gap_open)
{
    gap_scheme scheme{};
    EXPECT_EQ(scheme.get_gap_open_score(), 0);
    scheme.get_gap_open_score() = -2;
    EXPECT_EQ(scheme.get_gap_open_score(), -2);
}

TEST(gap_scheme, set_linear)
{
    gap_scheme scheme{gap_score{-2}};
    EXPECT_EQ(scheme.get_gap_score(), -2);
    EXPECT_EQ(scheme.get_gap_open_score(), 0);

    scheme.set_linear(gap_score{-3});
    EXPECT_EQ(scheme.get_gap_score(), -3);
    EXPECT_EQ(scheme.get_gap_open_score(), 0);
}

TEST(gap_scheme, set_affine)
{
    gap_scheme scheme{gap_score{-2}, gap_open_score{-4}};
    EXPECT_EQ(scheme.get_gap_score(), -2);
    EXPECT_EQ(scheme.get_gap_open_score(), -4);

    scheme.set_affine(gap_score{-3}, gap_open_score{-6});
    EXPECT_EQ(scheme.get_gap_score(), -3);
    EXPECT_EQ(scheme.get_gap_open_score(), -6);
}

TEST(gap_scheme, score)
{
    gap_scheme scheme{gap_score{-2}};
    EXPECT_EQ(scheme.score(0), 0);
    EXPECT_EQ(scheme.score(2), -4);
    EXPECT_EQ(scheme.score(5), -10);

    scheme.set_affine(gap_score{-3}, gap_open_score{-6});
    EXPECT_EQ(scheme.score(0), 0);
    EXPECT_EQ(scheme.score(2), -12);
    EXPECT_EQ(scheme.score(5), -21);
}

TEST(gap_scheme, serialisation)
{
    gap_scheme scheme1;

    scheme1.set_linear(gap_score{-3});
    test::do_serialisation(scheme1);

    scheme1.set_affine(gap_score{-3}, gap_open_score{-6});
    test::do_serialisation(scheme1);
}
