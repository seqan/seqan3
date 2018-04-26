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

#include <gtest/gtest.h>

#include <random>

#include <seqan3/std/concept/core_language.hpp>

#include "auxiliary.hpp"

using namespace seqan3;

TEST(same_concept, basic)
{
    EXPECT_TRUE((same_concept<int, int, int>));
    EXPECT_TRUE((!same_concept<int, char, int>));
}

TEST(derived_from_conept, basic)
{
    EXPECT_TRUE((derived_from_conept<type_b,
                                     type_a>));
    EXPECT_TRUE((!derived_from_conept<type_a,
                                      type_b>));
}

TEST(implicitly_convertible_to_concept, basic)
{
    EXPECT_TRUE((implicitly_convertible_to_concept<type_b,
                                                   type_c>));
    EXPECT_TRUE((!implicitly_convertible_to_concept<type_c,
                                                    type_b>));
    EXPECT_TRUE((!implicitly_convertible_to_concept<type_a,
                                                    type_c>));
}

TEST(explicitly_convertible_to_concept, basic)
{
    EXPECT_TRUE((explicitly_convertible_to_concept<type_b,
                                                   type_c>));
    EXPECT_TRUE((!explicitly_convertible_to_concept<type_c,
                                                    type_b>));
    EXPECT_TRUE((explicitly_convertible_to_concept<type_a,
                                                   type_c>));
}

TEST(convertible_to_concept, basic)
{
    EXPECT_TRUE((convertible_to_concept<type_b,
                                        type_c>));
    EXPECT_TRUE((!convertible_to_concept<type_c,
                                         type_b>));
    EXPECT_TRUE((!convertible_to_concept<type_a,
                                         type_c>));
}

TEST(common_reference_concept, basic)
{
    EXPECT_TRUE((common_reference_concept<int32_t, int16_t, int8_t>));
    EXPECT_TRUE((!common_reference_concept<int32_t, int16_t, type_c>));
}

TEST(common_concept, basic)
{
    EXPECT_TRUE((common_concept<type_a,
                                type_b>));
    EXPECT_TRUE((!common_reference_concept<type_a,
                                           type_c>));
}

TEST(integral_concept, basic)
{
    EXPECT_TRUE((integral_concept<int>));
    EXPECT_TRUE((!integral_concept<float>));
}

TEST(signed_integral_concept, basic)
{
    EXPECT_TRUE((signed_integral_concept<int>));
    EXPECT_TRUE((!signed_integral_concept<unsigned>));
}

TEST(unsigned_integral_concept, basic)
{
    EXPECT_TRUE((!unsigned_integral_concept<int>));
    EXPECT_TRUE((unsigned_integral_concept<unsigned>));
}

TEST(assignable_concept, basic)
{
    EXPECT_TRUE((assignable_concept<type_a&,
                                    type_a const &>));
    EXPECT_TRUE((assignable_concept<type_c&,
                                    type_b const &>));
    EXPECT_TRUE((!assignable_concept<type_a&,
                                     type_c&>));
}

TEST(swappable_concept, basic)
{
    EXPECT_TRUE((swappable_concept<type_a>));
    EXPECT_TRUE((swappable_concept<type_b>));
}

TEST(swappable_with_concept, basic)
{
    EXPECT_TRUE((swappable_with_concept<type_a&,
                                        type_a&>));
    EXPECT_TRUE((!swappable_with_concept<type_b,
                                         type_c>));
}
