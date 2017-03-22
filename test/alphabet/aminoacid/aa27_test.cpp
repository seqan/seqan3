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

#include <gtest/gtest.h>
#include <seqan3/alphabet/aminoacid/aa27.hpp>

using namespace seqan3;

TEST(alphabet_aminoacid_aa27_test, default_ctr)
{
    constexpr aa27 amino{};

    EXPECT_EQ(amino, aa27::A);
    EXPECT_EQ(amino.value, aa27::c_type::A);
}

TEST(alphabet_aminoacid_aa27_test, copy_ctr)
{
    constexpr aa27 amino1{aa27::B};
    constexpr aa27 amino2{amino1};

    EXPECT_EQ(amino2, aa27::B);
    EXPECT_EQ(amino2.value, aa27::c_type::B);
}

TEST(alphabet_aminoacid_aa27_test, move_ctr)
{
    constexpr aa27 amino1{aa27::C};
    constexpr aa27 amino2{std::move(amino1)};

    EXPECT_EQ(amino2, aa27::C);
    EXPECT_EQ(amino2.value, aa27::c_type::C);
}

TEST(alphabet_aminoacid_aa27_test, copy_asn)
{
    constexpr aa27 amino1{aa27::D};
    aa27 amino2{aa27::E};

    amino2 = amino1;

    EXPECT_EQ(amino2, aa27::D);
    EXPECT_EQ(amino2.value, aa27::c_type::D);
}

TEST(alphabet_aminoacid_aa27_test, move_asn)
{
    constexpr aa27 amino1{aa27::F};
    aa27 amino2{aa27::G};

    amino2 = std::move(amino1);

    EXPECT_EQ(amino2, aa27::F);
    EXPECT_EQ(amino2.value, aa27::c_type::F);
}

TEST(alphabet_aminoacid_aa27_test, operator_char)
{
    constexpr aa27 amino{aa27::G};
    constexpr char c = static_cast<char>(amino);

    EXPECT_EQ(c, 'G');
}

TEST(alphabet_aminoacid_aa27_test, to_char)
{
    constexpr aa27 amino_term{aa27::TERMINATOR};
    constexpr char c_term = amino_term.to_char();

    EXPECT_EQ(c_term, '*');

    for (auto i = 'A'; i <= 'Z'; ++i)
    {
        uint8_t integral = i - 'A';
        aa27 amino{aa27::c_type{integral}};
        char c = amino.to_char();
        EXPECT_EQ(c, i);
    }
}

TEST(alphabet_aminoacid_aa27_test, from_char)
{
    aa27 amino;
    amino = amino.from_char('*');

    EXPECT_EQ(amino, aa27::TERMINATOR);
    EXPECT_EQ(amino.value, aa27::c_type::TERMINATOR);

    for (auto i = 'A'; i <= 'Z'; ++i)
    {
        uint8_t integral = i - 'A';
        EXPECT_EQ(aa27{}.from_char(i), aa27::c_type{integral});
    }

    for (auto i = 'a'; i <= 'z'; ++i)
    {
        uint8_t integral = i - 'a';
        EXPECT_EQ(aa27{}.from_char(i), aa27::c_type{integral});
    }
}

TEST(alphabet_aminoacid_aa27_test, to_integral)
{
    constexpr  aa27 amino_term{aa27::TERMINATOR};
    constexpr  uint8_t integral_term = amino_term.to_integral();

    EXPECT_EQ(integral_term, 26);

    for (auto i = 'A'; i <= 'Z'; ++i)
    {
        uint8_t integral = i - 'A';
        aa27 amino{aa27::c_type{integral}};
        uint8_t integral_res = amino.to_integral();
        EXPECT_EQ(integral_res, integral);
    }

    for (auto i = 'a'; i <= 'z'; ++i)
    {
        uint8_t integral = i - 'a';
        aa27 amino{aa27::c_type{integral}};
        uint8_t integral_res = amino.to_integral();
        EXPECT_EQ(integral_res, integral);
    }
}

TEST(alphabet_aminoacid_aa27_test, from_integral)
{
    aa27 amino_term;
    amino_term = amino_term.from_integral(26);

    EXPECT_EQ(amino_term, aa27::TERMINATOR);
    EXPECT_EQ(amino_term.value, aa27::c_type::TERMINATOR);

    for (auto i = 'A'; i <= 'Z'; ++i)
    {
        uint8_t integral = i - 'A';
        aa27 amino;
        amino = amino.from_integral(integral);
        EXPECT_EQ(amino, aa27::c_type{integral});
    }
}

TEST(alphabet_aminoacid_aa27_test, logical_operations)
{
    EXPECT_TRUE(aa27{aa27::A} == aa27{aa27::A});
    EXPECT_TRUE(aa27{aa27::A} != aa27{aa27::T});
    EXPECT_TRUE(aa27{aa27::A} < aa27{aa27::C});
    EXPECT_TRUE(aa27{aa27::G} > aa27{aa27::C});
    EXPECT_TRUE(aa27{aa27::A} <= aa27{aa27::G});
    EXPECT_TRUE(aa27{aa27::G} >= aa27{aa27::A});
}
