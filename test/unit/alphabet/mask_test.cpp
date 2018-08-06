// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//

/*!\file
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \brief Test files for the mask composition.
 */


#include <gtest/gtest.h>

#include <seqan3/alphabet/concept_pre.hpp>
#include <seqan3/alphabet/mask/all.hpp>

using namespace seqan3;

TEST(assignment, assign_rank)
{
    // l-value
    mask lmask;
    EXPECT_EQ(lmask.assign_rank(1), mask::MASKED);
    EXPECT_TRUE(lmask.to_rank());
    EXPECT_EQ(lmask.assign_rank(0), mask::UNMASKED);
    EXPECT_FALSE(lmask.to_rank());
    EXPECT_EQ(lmask.assign_rank(true), mask::MASKED);
    EXPECT_EQ(lmask.assign_rank(false), mask::UNMASKED);
    EXPECT_EQ(lmask.assign_rank(mask::MASKED), mask::MASKED);
    EXPECT_TRUE(lmask.to_rank());
    EXPECT_EQ(lmask.assign_rank(mask::UNMASKED), mask::UNMASKED);
    EXPECT_FALSE(lmask.to_rank());

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
    EXPECT_EQ(std::move(rmask).assign_rank(mask::MASKED), mask::MASKED);
    EXPECT_TRUE(std::move(rmask).to_rank());
    EXPECT_EQ(std::move(rmask).assign_rank(mask::UNMASKED), mask::UNMASKED);
    EXPECT_FALSE(std::move(rmask).to_rank());

    // const r-value
    mask const crmask{lmask};
    EXPECT_EQ(std::move(crmask).to_rank(), lmask.to_rank());
    EXPECT_TRUE((std::is_same_v<decltype(std::move(crmask)), mask const &&>));
}

// ------------------------------------------------------------------
// comparators
// ------------------------------------------------------------------
TEST(comparators, compare)
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
