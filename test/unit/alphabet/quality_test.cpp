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
// ============================================================================

#include <numeric>
#include <sstream>
#include <tuple>
#include <utility>

#include <gtest/gtest.h>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/quality/all.hpp>

using namespace seqan3;

template <typename T>
class quality : public ::testing::Test
{};

// add all alphabets from the quality sub module here
using quality_types     = ::testing::Types<phred42, phred63, phred68legacy>;
using quality_types2    =       meta::list<phred42, phred63, phred68legacy>;

TYPED_TEST_CASE(quality, quality_types);

// more elaborate test of assign_char and to_char, basic test is in alphabet_test.cpp
TYPED_TEST(quality, conversion_char)
{
    using c_t = typename TypeParam::char_type;
    for (c_t i = std::numeric_limits<c_t>::lowest(); i < std::numeric_limits<c_t>::max(); ++i)
    {
        TypeParam v;
        v.assign_char(i);

        if (i < TypeParam::offset_char)                                     // too small, map to valid smallest
            EXPECT_EQ(v.to_char(), TypeParam::offset_char);
        else if (i >= TypeParam::offset_char + TypeParam::value_size)       // too big, map to valid biggest
            EXPECT_EQ(v.to_char(), TypeParam::offset_char + TypeParam::value_size - 1);
        else                                                                // valid range, map to identity
            EXPECT_EQ(v.to_char(), i);
    }
}

// test assign_phred and to_phred
TYPED_TEST(quality, conversion_phred)
{
    using p_t = typename TypeParam::phred_type;
    for (p_t i = std::numeric_limits<p_t>::lowest(); i < std::numeric_limits<p_t>::max(); ++i)
    {
        TypeParam v;
        v.assign_phred(i);

        if (i < TypeParam::offset_phred)                                     // too small, map to valid smallest
            EXPECT_EQ(v.to_phred(), TypeParam::offset_phred);
        else if (i >= TypeParam::offset_phred + TypeParam::value_size)       // too big, map to valid biggest
            EXPECT_EQ(v.to_phred(), TypeParam::offset_phred + TypeParam::value_size - 1);
        else                                                                // valid range, map to identity
            EXPECT_EQ(v.to_phred(), i);
    }
}

// test user-defined constructor
TYPED_TEST(quality, construction_by_phred)
{
    TypeParam v{0};
    EXPECT_EQ(v.to_phred(), 0);
    EXPECT_EQ(v.to_rank(),  -TypeParam::offset_phred);

    TypeParam v2{23};
    EXPECT_EQ(v2.to_phred(), 23);
    EXPECT_EQ(v2.to_rank(),  23 - TypeParam::offset_phred);
}

// ------------------------------------------------------------------
// conversion
// ------------------------------------------------------------------

// test provision of data type `phred_type` and phred converter.
TYPED_TEST(quality, quality_concept)
{
    EXPECT_TRUE(quality_concept<TypeParam>);
}

TYPED_TEST(quality, explicit_conversion)
{
    meta::for_each(quality_types2{}, [&] (auto && qual) constexpr
    {
        using out_type = std::decay_t<decltype(qual)>;
        EXPECT_EQ(static_cast<out_type>(TypeParam{ 0}), out_type{ 0});
        EXPECT_EQ(static_cast<out_type>(TypeParam{ 5}), out_type{ 5});
        EXPECT_EQ(static_cast<out_type>(TypeParam{15}), out_type{15});
        EXPECT_EQ(static_cast<out_type>(TypeParam{20}), out_type{20});
        EXPECT_EQ(static_cast<out_type>(TypeParam{40}), out_type{40});
    });
}
