// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
// Author: Sara Hetzel <sara.hetzel AT fu-berlin.de>
// ============================================================================

#include <array>

#include <gtest/gtest.h>
#include <seqan3/alphabet/nucleotide/nucl16.hpp>

using namespace seqan3;

constexpr char capital_char[] = {'A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R', 'S', 'T', 'U', 'V', 'W', 'Y'};
constexpr char lower_char[] = {'a', 'b', 'c', 'd', 'g', 'h', 'k', 'm', 'n', 'r', 's', 't', 'u', 'v', 'w', 'y'};
constexpr nucl16 all_nucl16[] = {nucl16::A, nucl16::B, nucl16::C, nucl16::D, nucl16::G, nucl16::H, nucl16::K,
                                 nucl16::M, nucl16::N, nucl16::R, nucl16::S, nucl16::T, nucl16::U, nucl16::V,
                                 nucl16::W, nucl16::Y};

TEST(alphabet_nucleotides_nucl16_test, default_ctr)
{
    constexpr nucl16 nucleotide{};

    EXPECT_EQ(nucleotide, nucl16::A);
}

TEST(alphabet_nucleotides_nucl16_test, copy_ctr)
{
    constexpr nucl16 nucleotide1{nucl16::B};
    constexpr nucl16 nucleotide2{nucleotide1};

    EXPECT_EQ(nucleotide2, nucl16::B);
}

TEST(alphabet_nucleotides_nucl16_test, move_ctr)
{
    constexpr nucl16 nucleotide1{nucl16::C};
    constexpr nucl16 nucleotide2{std::move(nucleotide1)};

    EXPECT_EQ(nucleotide2, nucl16::C);
}

TEST(alphabet_nucleotides_nucl16_test, copy_asn)
{
    constexpr nucl16 nucleotide1{nucl16::D};
    nucl16 nucleotide2{nucl16::H};

    nucleotide2 = nucleotide1;

    EXPECT_EQ(nucleotide2, nucl16::D);
}

TEST(alphabet_nucleotides_nucl16_test, move_asn)
{
    constexpr nucl16 nucleotide1{nucl16::K};
    nucl16 nucleotide2{nucl16::G};

    nucleotide2 = std::move(nucleotide1);

    EXPECT_EQ(nucleotide2, nucl16::K);
}

TEST(alphabet_nucleotides_nucl16_test, operator_char)
{
    constexpr nucl16 nucleotide{nucl16::G};
    constexpr char c = static_cast<char>(nucleotide);

    EXPECT_EQ(c, 'G');
}

TEST(alphabet_nucleotides_nucl16_test, to_char)
{
    uint8_t rank = 0;
    for (auto i : capital_char)
    {
        nucl16 nucleotide{all_nucl16[rank++]};
        char c = nucleotide.to_char();
        EXPECT_EQ(c, i);
    }
}

TEST(alphabet_nucleotides_nucl16_test, from_char)
{
    uint8_t rank = 0;
    for (auto i : capital_char)
    {
        nucl16 nucleotide{all_nucl16[rank++]};
        EXPECT_EQ(nucl16{}.from_char(i), nucleotide);
    }

    rank = 0;
    for (auto i : lower_char)
    {
        nucl16 nucleotide{all_nucl16[rank++]};
        EXPECT_EQ(nucl16{}.from_char(i), nucleotide);
    }
}

TEST(alphabet_nucleotides_nucl16_test, to_integral)
{
    uint8_t rank = 0;
    for (auto i : capital_char)
    {
        uint8_t integral = rank++;
        nucl16 nucleotide{all_nucl16[integral]};
        uint8_t integral_res = nucleotide.to_integral();
        EXPECT_EQ(integral_res, integral);
    }

    rank = 0;
    for (auto i : lower_char)
    {
        uint8_t integral = rank++;
        nucl16 nucleotide{all_nucl16[integral]};
        uint8_t integral_res = nucleotide.to_integral();
        EXPECT_EQ(integral_res, integral);
    }
}

TEST(alphabet_nucleotides_nucl16_test, from_integral)
{
    uint8_t rank = 0;
    for (auto i : capital_char)
    {
        uint8_t integral = rank++;
        nucl16 nucleotide;
        nucleotide = nucleotide.from_integral(integral);
        EXPECT_EQ(nucleotide, all_nucl16[integral]);
    }
}

TEST(alphabet_nucleotides_nucl16_test, logical_operations)
{
    EXPECT_TRUE(nucl16{nucl16::A} == nucl16{nucl16::A});
    EXPECT_TRUE(nucl16{nucl16::A} != nucl16{nucl16::T});
    EXPECT_TRUE(nucl16{nucl16::A} < nucl16{nucl16::C});
    EXPECT_TRUE(nucl16{nucl16::G} > nucl16{nucl16::C});
    EXPECT_TRUE(nucl16{nucl16::A} <= nucl16{nucl16::G});
    EXPECT_TRUE(nucl16{nucl16::G} >= nucl16{nucl16::A});
}
