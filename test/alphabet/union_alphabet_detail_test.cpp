// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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
// ==========================================================================

#include <gtest/gtest.h>

#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/union_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

using namespace seqan3;

TEST(union_alphabet_detail_test, alphabet_prefix_sum_sizes)
{
    constexpr auto sizes0 = detail::alphabet_prefix_sum_sizes<>();
    EXPECT_EQ(sizes0.size(), 1);
    EXPECT_EQ(sizes0[0], 0);

    constexpr auto sizes1 = detail::alphabet_prefix_sum_sizes<dna4>();
    EXPECT_EQ(sizes1.size(), 2);
    EXPECT_EQ(sizes1[0], 0);
    EXPECT_EQ(sizes1[1], 4);

    constexpr auto sizes2 = detail::alphabet_prefix_sum_sizes<dna4, gap>();
    EXPECT_EQ(sizes2.size(), 3);
    EXPECT_EQ(sizes2[0], 0);
    EXPECT_EQ(sizes2[1], 4);
    EXPECT_EQ(sizes2[2], 5);

    constexpr auto sizes3 = detail::alphabet_prefix_sum_sizes<dna4, gap, dna5>();
    EXPECT_EQ(sizes3.size(), 4);
    EXPECT_EQ(sizes3[0], 0);
    EXPECT_EQ(sizes3[1], 4);
    EXPECT_EQ(sizes3[2], 5);
    EXPECT_EQ(sizes3[3], 10);

    constexpr auto sizes4 = detail::alphabet_prefix_sum_sizes<dna5, gap, dna4>();
    EXPECT_EQ(sizes4.size(), 4);
    EXPECT_EQ(sizes4[0], 0);
    EXPECT_EQ(sizes4[1], 5);
    EXPECT_EQ(sizes4[2], 6);
    EXPECT_EQ(sizes4[3], 10);
}

TEST(union_alphabet_detail_test, union_alphabet_value_to_char_tableI)
{
    constexpr auto table1 = detail::union_alphabet::value_to_char_table_I<5, char>(dna4{});
    EXPECT_EQ(table1.size(), 5);
    EXPECT_EQ(table1[0], 'A');
    EXPECT_EQ(table1[1], 'C');
    EXPECT_EQ(table1[2], 'G');
    EXPECT_EQ(table1[3], 'T');
    EXPECT_EQ(table1[4], '\0');

    constexpr auto table2 = detail::union_alphabet::value_to_char_table_I<5, char>(dna5{});
    EXPECT_EQ(table2.size(), 5);
    EXPECT_EQ(table2[0], 'A');
    EXPECT_EQ(table2[1], 'C');
    EXPECT_EQ(table2[2], 'G');
    EXPECT_EQ(table2[3], 'T');
    EXPECT_EQ(table2[4], 'N');

    constexpr auto table3 = detail::union_alphabet::value_to_char_table_I<5, char>(gap{});
    EXPECT_EQ(table3.size(), 5);
    EXPECT_EQ(table3[0], '-');
    EXPECT_EQ(table3[1], '\0');
    EXPECT_EQ(table3[2], '\0');
    EXPECT_EQ(table3[3], '\0');
    EXPECT_EQ(table3[4], '\0');
}

TEST(union_alphabet_detail_test, union_alphabet_value_to_char_table)
{
    constexpr auto value_to_char0 = detail::union_alphabet::value_to_char_table<char>();
    EXPECT_EQ(value_to_char0.size(), 0);

    constexpr auto value_to_char1 = detail::union_alphabet::value_to_char_table<char, dna4>();
    EXPECT_EQ(value_to_char1.size(), 4);
    EXPECT_EQ(value_to_char1[0], 'A');
    EXPECT_EQ(value_to_char1[1], 'C');
    EXPECT_EQ(value_to_char1[2], 'G');
    EXPECT_EQ(value_to_char1[3], 'T');

    constexpr auto value_to_char2 = detail::union_alphabet::value_to_char_table<char, dna4, gap>();
    EXPECT_EQ(value_to_char2.size(), 5);
    EXPECT_EQ(value_to_char2[0], 'A');
    EXPECT_EQ(value_to_char2[1], 'C');
    EXPECT_EQ(value_to_char2[2], 'G');
    EXPECT_EQ(value_to_char2[3], 'T');
    EXPECT_EQ(value_to_char2[4], '-');

    constexpr auto value_to_char3 = detail::union_alphabet::value_to_char_table<char, dna4, gap, dna5>();
    EXPECT_EQ(value_to_char3.size(), 10);
    EXPECT_EQ(value_to_char3[0], 'A');
    EXPECT_EQ(value_to_char3[1], 'C');
    EXPECT_EQ(value_to_char3[2], 'G');
    EXPECT_EQ(value_to_char3[3], 'T');
    EXPECT_EQ(value_to_char3[4], '-');
    EXPECT_EQ(value_to_char3[5], 'A');
    EXPECT_EQ(value_to_char3[6], 'C');
    EXPECT_EQ(value_to_char3[7], 'G');
    EXPECT_EQ(value_to_char3[8], 'T');
    EXPECT_EQ(value_to_char3[9], 'N');

    constexpr auto value_to_char4 = detail::union_alphabet::value_to_char_table<char, dna5, gap, dna4>();
    EXPECT_EQ(value_to_char4.size(), 10);
    EXPECT_EQ(value_to_char4[0], 'A');
    EXPECT_EQ(value_to_char4[1], 'C');
    EXPECT_EQ(value_to_char4[2], 'G');
    EXPECT_EQ(value_to_char4[3], 'T');
    EXPECT_EQ(value_to_char4[4], 'N');
    EXPECT_EQ(value_to_char4[5], '-');
    EXPECT_EQ(value_to_char4[6], 'A');
    EXPECT_EQ(value_to_char4[7], 'C');
    EXPECT_EQ(value_to_char4[8], 'G');
    EXPECT_EQ(value_to_char4[9], 'T');
}

TEST(union_alphabet_detail_test, union_alphabet_char_to_value_table)
{
    constexpr auto char_to_value0 = detail::union_alphabet::char_to_value_table<char>();
    EXPECT_EQ(char_to_value0.size(), 256);
    for(auto value: char_to_value0)
        EXPECT_EQ(value, 0);

    constexpr auto char_to_value1 = detail::union_alphabet::char_to_value_table<char, dna4>();
    EXPECT_EQ(char_to_value1.size(), 256);
    EXPECT_EQ(char_to_value1['A'], 0);
    EXPECT_EQ(char_to_value1['C'], 1);
    EXPECT_EQ(char_to_value1['G'], 2);
    EXPECT_EQ(char_to_value1['T'], 3);

    constexpr auto char_to_value2 = detail::union_alphabet::char_to_value_table<char, dna4, gap>();
    EXPECT_EQ(char_to_value2.size(), 256);
    EXPECT_EQ(char_to_value2['A'], 0);
    EXPECT_EQ(char_to_value2['C'], 1);
    EXPECT_EQ(char_to_value2['G'], 2);
    EXPECT_EQ(char_to_value2['T'], 3);
    EXPECT_EQ(char_to_value2['-'], 4);

    constexpr auto char_to_value3 = detail::union_alphabet::char_to_value_table<char, dna4, gap, dna5>();
    EXPECT_EQ(char_to_value3.size(), 256);
    EXPECT_EQ(char_to_value3['A'], 0);
    EXPECT_EQ(char_to_value3['C'], 1);
    EXPECT_EQ(char_to_value3['G'], 2);
    EXPECT_EQ(char_to_value3['T'], 3);
    EXPECT_EQ(char_to_value3['-'], 4);
    EXPECT_EQ(char_to_value3['A'], 0);
    EXPECT_EQ(char_to_value3['C'], 1);
    EXPECT_EQ(char_to_value3['G'], 2);
    EXPECT_EQ(char_to_value3['T'], 3);
    EXPECT_EQ(char_to_value3['N'], 9);

    constexpr auto char_to_value4 = detail::union_alphabet::char_to_value_table<char, dna5, gap, dna4>();
    EXPECT_EQ(char_to_value4.size(), 256);
    EXPECT_EQ(char_to_value4['A'], 0);
    EXPECT_EQ(char_to_value4['C'], 1);
    EXPECT_EQ(char_to_value4['G'], 2);
    EXPECT_EQ(char_to_value4['T'], 3);
    EXPECT_EQ(char_to_value4['N'], 4);
    EXPECT_EQ(char_to_value4['-'], 5);
    EXPECT_EQ(char_to_value4['A'], 0);
    EXPECT_EQ(char_to_value4['C'], 1);
    EXPECT_EQ(char_to_value4['G'], 2);
    EXPECT_EQ(char_to_value4['T'], 3);
}
