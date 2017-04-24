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

constexpr aa27 all_aa27[] = {aa27::A, aa27::B, aa27::C, aa27::D, aa27::E, aa27::F, aa27::G, aa27::H, aa27::I, aa27::J,
                             aa27::K, aa27::L, aa27::M, aa27::N, aa27::O, aa27::P, aa27::Q, aa27::R, aa27::S, aa27::T,
                             aa27::U, aa27::V, aa27::W, aa27::X, aa27::Y, aa27::Z};

TEST(alphabet_aminoacid_aa27_test, default_ctr)
{
    constexpr aa27 amino{};

    EXPECT_EQ(amino, aa27::A);
}

TEST(alphabet_aminoacid_aa27_test, copy_ctr)
{
    constexpr aa27 amino1{aa27::B};
    constexpr aa27 amino2{amino1};

    EXPECT_EQ(amino2, aa27::B);
}

TEST(alphabet_aminoacid_aa27_test, move_ctr)
{
    constexpr aa27 amino1{aa27::C};
    constexpr aa27 amino2{std::move(amino1)};

    EXPECT_EQ(amino2, aa27::C);
}

TEST(alphabet_aminoacid_aa27_test, copy_asn)
{
    constexpr aa27 amino1{aa27::D};
    aa27 amino2{aa27::E};

    amino2 = amino1;

    EXPECT_EQ(amino2, aa27::D);
}

TEST(alphabet_aminoacid_aa27_test, move_asn)
{
    constexpr aa27 amino1{aa27::F};
    aa27 amino2{aa27::G};

    amino2 = std::move(amino1);

    EXPECT_EQ(amino2, aa27::F);
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
        uint8_t rank = i - 'A';
        aa27 amino{all_aa27[rank]};
        char c = amino.to_char();
        EXPECT_EQ(c, i);
    }
}

TEST(alphabet_aminoacid_aa27_test, assign_char)
{
    aa27 amino;
    amino = amino.assign_char('*');

    EXPECT_EQ(amino, aa27::TERMINATOR);

    for (auto i = 'A'; i <= 'Z'; ++i)
    {
        uint8_t rank = i - 'A';
        EXPECT_EQ(aa27{}.assign_char(i), all_aa27[rank]);
    }

    for (auto i = 'a'; i <= 'z'; ++i)
    {
        uint8_t rank = i - 'a';
        EXPECT_EQ(aa27{}.assign_char(i), all_aa27[rank]);
    }
}

TEST(alphabet_aminoacid_aa27_test, to_rank)
{
    constexpr  aa27 amino_term{aa27::TERMINATOR};
    constexpr  uint8_t rank_term = amino_term.to_rank();

    EXPECT_EQ(rank_term, 26);

    for (auto i = 'A'; i <= 'Z'; ++i)
    {
        uint8_t rank = i - 'A';
        aa27 amino{all_aa27[rank]};
        uint8_t rank_res = amino.to_rank();
        EXPECT_EQ(rank_res, rank);
    }

    for (auto i = 'a'; i <= 'z'; ++i)
    {
        uint8_t rank = i - 'a';
        aa27 amino{all_aa27[rank]};
        uint8_t rank_res = amino.to_rank();
        EXPECT_EQ(rank_res, rank);
    }
}

TEST(alphabet_aminoacid_aa27_test, assign_rank)
{
    aa27 amino_term;
    amino_term = amino_term.assign_rank(26);

    EXPECT_EQ(amino_term, aa27::TERMINATOR);

    for (auto i = 'A'; i <= 'Z'; ++i)
    {
        uint8_t rank = i - 'A';
        aa27 amino;
        amino = amino.assign_rank(rank);
        EXPECT_EQ(amino, all_aa27[rank]);
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
