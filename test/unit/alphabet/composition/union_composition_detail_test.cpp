// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// ==========================================================================

#include <type_traits>

#include <gtest/gtest.h>

#include <seqan3/alphabet/composition/union_composition.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

using namespace seqan3;

template <typename ...alphabet_types>
class detail_union_composition : public union_composition<alphabet_types...>
{
public:
    using union_composition<alphabet_types...>::partial_sum_sizes;
    using union_composition<alphabet_types...>::rank_to_char;
    using union_composition<alphabet_types...>::char_to_rank;
};

TEST(union_composition_detail_test, partial_sum_sizes)
{
    constexpr std::array partial_sum2 = detail_union_composition<dna4, gap>::partial_sum_sizes;
    EXPECT_EQ(partial_sum2.size(), 3u);
    EXPECT_EQ(partial_sum2[0], 0);
    EXPECT_EQ(partial_sum2[1], 4);
    EXPECT_EQ(partial_sum2[2], 5);

    constexpr std::array partial_sum3 = detail_union_composition<dna4, gap, dna5>::partial_sum_sizes;
    EXPECT_EQ(partial_sum3.size(), 4u);
    EXPECT_EQ(partial_sum3[0], 0);
    EXPECT_EQ(partial_sum3[1], 4);
    EXPECT_EQ(partial_sum3[2], 5);
    EXPECT_EQ(partial_sum3[3], 10);

    constexpr std::array partial_sum4 = detail_union_composition<dna5, gap, dna4>::partial_sum_sizes;
    EXPECT_EQ(partial_sum4.size(), 4u);
    EXPECT_EQ(partial_sum4[0], 0);
    EXPECT_EQ(partial_sum4[1], 5);
    EXPECT_EQ(partial_sum4[2], 6);
    EXPECT_EQ(partial_sum4[3], 10);
}

TEST(union_composition_detail_test, rank_to_char)
{
    constexpr std::array rank_to_char2 = detail_union_composition<dna4, gap>::rank_to_char;
    EXPECT_EQ(rank_to_char2.size(), 5u);
    EXPECT_EQ(rank_to_char2[0], 'A');
    EXPECT_EQ(rank_to_char2[1], 'C');
    EXPECT_EQ(rank_to_char2[2], 'G');
    EXPECT_EQ(rank_to_char2[3], 'T');
    EXPECT_EQ(rank_to_char2[4], '-');

    constexpr std::array rank_to_char3 = detail_union_composition<dna4, gap, dna5>::rank_to_char;
    EXPECT_EQ(rank_to_char3.size(), 10u);
    EXPECT_EQ(rank_to_char3[0], 'A');
    EXPECT_EQ(rank_to_char3[1], 'C');
    EXPECT_EQ(rank_to_char3[2], 'G');
    EXPECT_EQ(rank_to_char3[3], 'T');
    EXPECT_EQ(rank_to_char3[4], '-');
    EXPECT_EQ(rank_to_char3[5], 'A');
    EXPECT_EQ(rank_to_char3[6], 'C');
    EXPECT_EQ(rank_to_char3[7], 'G');
    EXPECT_EQ(rank_to_char3[8], 'N');
    EXPECT_EQ(rank_to_char3[9], 'T');

    constexpr std::array rank_to_char4 = detail_union_composition<dna5, gap, dna4>::rank_to_char;
    EXPECT_EQ(rank_to_char4.size(), 10u);
    EXPECT_EQ(rank_to_char4[0], 'A');
    EXPECT_EQ(rank_to_char4[1], 'C');
    EXPECT_EQ(rank_to_char4[2], 'G');
    EXPECT_EQ(rank_to_char4[3], 'N');
    EXPECT_EQ(rank_to_char4[4], 'T');
    EXPECT_EQ(rank_to_char4[5], '-');
    EXPECT_EQ(rank_to_char4[6], 'A');
    EXPECT_EQ(rank_to_char4[7], 'C');
    EXPECT_EQ(rank_to_char4[8], 'G');
    EXPECT_EQ(rank_to_char4[9], 'T');
}

TEST(union_composition_detail_test, char_to_rank)
{
    constexpr std::array char_to_rank2 = detail_union_composition<dna4, gap>::char_to_rank;
    EXPECT_EQ(char_to_rank2.size(), 256u);
    EXPECT_EQ(char_to_rank2['A'], 0);
    EXPECT_EQ(char_to_rank2['C'], 1);
    EXPECT_EQ(char_to_rank2['G'], 2);
    EXPECT_EQ(char_to_rank2['T'], 3);
    EXPECT_EQ(char_to_rank2['-'], 4);

    constexpr std::array char_to_rank3 = detail_union_composition<dna4, gap, dna5>::char_to_rank;
    EXPECT_EQ(char_to_rank3.size(), 256u);
    EXPECT_EQ(char_to_rank3['A'], 0);
    EXPECT_EQ(char_to_rank3['C'], 1);
    EXPECT_EQ(char_to_rank3['G'], 2);
    EXPECT_EQ(char_to_rank3['T'], 3);
    EXPECT_EQ(char_to_rank3['-'], 4);
    EXPECT_EQ(char_to_rank3['A'], 0);
    EXPECT_EQ(char_to_rank3['C'], 1);
    EXPECT_EQ(char_to_rank3['G'], 2);
    EXPECT_EQ(char_to_rank3['N'], 8);
    EXPECT_EQ(char_to_rank3['T'], 3);

    constexpr std::array char_to_rank4 = detail_union_composition<dna5, gap, dna4>::char_to_rank;
    EXPECT_EQ(char_to_rank4.size(), 256u);
    EXPECT_EQ(char_to_rank4['A'], 0);
    EXPECT_EQ(char_to_rank4['C'], 1);
    EXPECT_EQ(char_to_rank4['G'], 2);
    EXPECT_EQ(char_to_rank4['N'], 3);
    EXPECT_EQ(char_to_rank4['T'], 4);
    EXPECT_EQ(char_to_rank4['-'], 5);
    EXPECT_EQ(char_to_rank4['A'], 0);
    EXPECT_EQ(char_to_rank4['C'], 1);
    EXPECT_EQ(char_to_rank4['G'], 2);
    EXPECT_EQ(char_to_rank4['T'], 4);
}
