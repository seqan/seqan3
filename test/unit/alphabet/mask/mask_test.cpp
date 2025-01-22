// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/mask/mask.hpp>

#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

INSTANTIATE_TYPED_TEST_SUITE_P(mask, semi_alphabet_test, seqan3::mask, );
INSTANTIATE_TYPED_TEST_SUITE_P(mask, semi_alphabet_constexpr, seqan3::mask, );

TEST(mask, assign_rank)
{
    // l-value
    seqan3::mask lmask;
    EXPECT_EQ(lmask.assign_rank(1), seqan3::mask::masked);
    EXPECT_TRUE(lmask.to_rank());
    EXPECT_EQ(lmask.assign_rank(0), seqan3::mask::unmasked);
    EXPECT_FALSE(lmask.to_rank());
    EXPECT_EQ(lmask.assign_rank(true), seqan3::mask::masked);
    EXPECT_EQ(lmask.assign_rank(false), seqan3::mask::unmasked);

    // const l-value
    lmask.assign_rank(1);
    seqan3::mask const clmask{lmask};
    EXPECT_TRUE(clmask.to_rank());

    // r-value
    seqan3::mask rmask{lmask};
    EXPECT_EQ(std::move(rmask).to_rank(), lmask.to_rank());
    EXPECT_TRUE((std::is_same_v<decltype(std::move(rmask)), seqan3::mask &&>));
    EXPECT_EQ(std::move(rmask).assign_rank(1), seqan3::mask::masked);
    EXPECT_TRUE(std::move(rmask).to_rank());
    EXPECT_EQ(std::move(rmask).assign_rank(0), seqan3::mask::unmasked);
    EXPECT_FALSE(std::move(rmask).to_rank());
    EXPECT_EQ(std::move(rmask).assign_rank(true), seqan3::mask::masked);
    EXPECT_EQ(std::move(rmask).assign_rank(false), seqan3::mask::unmasked);

    // const r-value
    seqan3::mask const crmask{lmask};
    EXPECT_EQ(std::move(crmask).to_rank(), lmask.to_rank());
    EXPECT_TRUE((std::is_same_v<decltype(std::move(crmask)), seqan3::mask const &&>));
}
