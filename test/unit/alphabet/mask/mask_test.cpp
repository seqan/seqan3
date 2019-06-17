// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/mask/mask.hpp>

using namespace seqan3;

TEST(mask, concept_check)
{
    EXPECT_TRUE(Semialphabet<mask>);
    EXPECT_TRUE(detail::ConstexprSemialphabet<mask>);
}

TEST(mask, assign_rank)
{
    // l-value
    mask lmask;
    EXPECT_EQ(lmask.assign_rank(1), mask::MASKED);
    EXPECT_TRUE(lmask.to_rank());
    EXPECT_EQ(lmask.assign_rank(0), mask::UNMASKED);
    EXPECT_FALSE(lmask.to_rank());
    EXPECT_EQ(lmask.assign_rank(true), mask::MASKED);
    EXPECT_EQ(lmask.assign_rank(false), mask::UNMASKED);

    // const l-value
    lmask.assign_rank(1);
    mask const clmask{lmask};
    EXPECT_TRUE(clmask.to_rank());

    // r-value
    mask rmask{lmask};
    EXPECT_EQ(std::move(rmask).to_rank(), lmask.to_rank());
    EXPECT_TRUE((std::is_same_v<decltype(std::move(rmask)), mask &&>));
    EXPECT_EQ(std::move(rmask).assign_rank(1), mask::MASKED);
    EXPECT_TRUE(std::move(rmask).to_rank());
    EXPECT_EQ(std::move(rmask).assign_rank(0), mask::UNMASKED);
    EXPECT_FALSE(std::move(rmask).to_rank());
    EXPECT_EQ(std::move(rmask).assign_rank(true), mask::MASKED);
    EXPECT_EQ(std::move(rmask).assign_rank(false), mask::UNMASKED);

    // const r-value
    mask const crmask{lmask};
    EXPECT_EQ(std::move(crmask).to_rank(), lmask.to_rank());
    EXPECT_TRUE((std::is_same_v<decltype(std::move(crmask)), mask const &&>));
}

// ------------------------------------------------------------------
// comparators
// ------------------------------------------------------------------
TEST(mask, compare)
{
    EXPECT_TRUE(mask::MASKED == mask::MASKED);
    EXPECT_TRUE(mask::MASKED != mask::UNMASKED);
    EXPECT_TRUE(mask::MASKED > mask::UNMASKED);
    EXPECT_TRUE(mask::UNMASKED < mask::MASKED);
    EXPECT_TRUE(mask::MASKED >= mask::UNMASKED);
    EXPECT_TRUE(mask::UNMASKED <= mask::MASKED);

    EXPECT_FALSE(mask::MASKED == mask::UNMASKED);
    EXPECT_FALSE(mask::MASKED != mask::MASKED);
    EXPECT_FALSE(mask::MASKED > mask::MASKED);
    EXPECT_FALSE(mask::UNMASKED < mask::UNMASKED);
    EXPECT_FALSE(mask::MASKED <= mask::UNMASKED);
    EXPECT_FALSE(mask::UNMASKED >= mask::MASKED);
}
