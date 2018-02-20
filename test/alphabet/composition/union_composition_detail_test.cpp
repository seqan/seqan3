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
    using union_composition<alphabet_types...>::sum_of_alphabet_sizes_v;
    using union_composition<alphabet_types...>::max_of_alphabet_sizes_v;
    using union_composition<alphabet_types...>::partial_sum_of_alphabet_sizes;
    using union_composition<alphabet_types...>::value_to_char_table_I;
    using union_composition<alphabet_types...>::value_to_char_table;
    using union_composition<alphabet_types...>::char_to_value_table;
};

TEST(union_composition_detail_test, sum_of_alphabet_sizes_v)
{
    using detail = detail_union_composition<dna4>;

    constexpr size_t sum0 = detail::sum_of_alphabet_sizes_v<>;
    constexpr size_t sum1 = detail::sum_of_alphabet_sizes_v<dna4>;
    constexpr size_t sum2 = detail::sum_of_alphabet_sizes_v<dna4, gap>;
    constexpr size_t sum3 = detail::sum_of_alphabet_sizes_v<dna4, gap, dna5>;
    constexpr size_t sum4 = detail::sum_of_alphabet_sizes_v<dna5, gap, dna4>;

    EXPECT_EQ(sum0, 0u);
    EXPECT_EQ(sum1, 4u);
    EXPECT_EQ(sum2, 5u);
    EXPECT_EQ(sum3, 10u);
    EXPECT_EQ(sum4, 10u);
}

TEST(union_composition_detail_test, max_of_alphabet_sizes_v)
{
    using detail = detail_union_composition<dna4>;

    constexpr size_t max0 = detail::max_of_alphabet_sizes_v<>;
    constexpr size_t max1 = detail::max_of_alphabet_sizes_v<dna4>;
    constexpr size_t max2 = detail::max_of_alphabet_sizes_v<dna4, gap>;
    constexpr size_t max3 = detail::max_of_alphabet_sizes_v<dna4, gap, dna5>;
    constexpr size_t max4 = detail::max_of_alphabet_sizes_v<dna5, gap, dna4>;

    EXPECT_EQ(max0, 0u);
    EXPECT_EQ(max1, 4u);
    EXPECT_EQ(max2, 4u);
    EXPECT_EQ(max3, 5u);
    EXPECT_EQ(max4, 5u);
}

TEST(union_composition_detail_test, partial_sum_of_alphabet_sizes)
{
    using detail = detail_union_composition<dna4>;

    constexpr std::array partial_sum0 = detail::partial_sum_of_alphabet_sizes<>();
    EXPECT_EQ(partial_sum0.size(), 1u);
    EXPECT_EQ(partial_sum0[0], 0);

    constexpr std::array partial_sum1 = detail::partial_sum_of_alphabet_sizes<dna4>();
    EXPECT_EQ(partial_sum1.size(), 2u);
    EXPECT_EQ(partial_sum1[0], 0);
    EXPECT_EQ(partial_sum1[1], 4);

    constexpr std::array partial_sum2 = detail::partial_sum_of_alphabet_sizes<dna4, gap>();
    EXPECT_EQ(partial_sum2.size(), 3u);
    EXPECT_EQ(partial_sum2[0], 0);
    EXPECT_EQ(partial_sum2[1], 4);
    EXPECT_EQ(partial_sum2[2], 5);

    constexpr std::array partial_sum3 = detail::partial_sum_of_alphabet_sizes<dna4, gap, dna5>();
    EXPECT_EQ(partial_sum3.size(), 4u);
    EXPECT_EQ(partial_sum3[0], 0);
    EXPECT_EQ(partial_sum3[1], 4);
    EXPECT_EQ(partial_sum3[2], 5);
    EXPECT_EQ(partial_sum3[3], 10);

    constexpr std::array partial_sum4 = detail::partial_sum_of_alphabet_sizes<dna5, gap, dna4>();
    EXPECT_EQ(partial_sum4.size(), 4u);
    EXPECT_EQ(partial_sum4[0], 0);
    EXPECT_EQ(partial_sum4[1], 5);
    EXPECT_EQ(partial_sum4[2], 6);
    EXPECT_EQ(partial_sum4[3], 10);
}

TEST(union_composition_detail_test, union_composition_value_to_char_table_I)
{
    using detail = detail_union_composition<dna4>;
    EXPECT_TRUE((std::is_same_v<detail::char_type, char>));

    constexpr std::array table1 = detail::value_to_char_table_I<5>(dna4{});
    EXPECT_EQ(table1.size(), 5u);
    EXPECT_EQ(table1[0], 'A');
    EXPECT_EQ(table1[1], 'C');
    EXPECT_EQ(table1[2], 'G');
    EXPECT_EQ(table1[3], 'T');
    EXPECT_EQ(table1[4], '\0');

    constexpr std::array table2 = detail::value_to_char_table_I<5>(dna5{});
    EXPECT_EQ(table2.size(), 5u);
    EXPECT_EQ(table2[0], 'A');
    EXPECT_EQ(table2[1], 'C');
    EXPECT_EQ(table2[2], 'G');
    EXPECT_EQ(table2[3], 'T');
    EXPECT_EQ(table2[4], 'N');

    constexpr std::array table3 = detail::value_to_char_table_I<5>(gap{});
    EXPECT_EQ(table3.size(), 5u);
    EXPECT_EQ(table3[0], '-');
    EXPECT_EQ(table3[1], '\0');
    EXPECT_EQ(table3[2], '\0');
    EXPECT_EQ(table3[3], '\0');
    EXPECT_EQ(table3[4], '\0');
}

TEST(union_composition_detail_test, union_composition_value_to_char_table)
{
    constexpr std::array value_to_char1 = detail_union_composition<dna4>::value_to_char_table();
    EXPECT_EQ(value_to_char1.size(), 4u);
    EXPECT_EQ(value_to_char1[0], 'A');
    EXPECT_EQ(value_to_char1[1], 'C');
    EXPECT_EQ(value_to_char1[2], 'G');
    EXPECT_EQ(value_to_char1[3], 'T');

    constexpr std::array value_to_char2 = detail_union_composition<dna4, gap>::value_to_char_table();
    EXPECT_EQ(value_to_char2.size(), 5u);
    EXPECT_EQ(value_to_char2[0], 'A');
    EXPECT_EQ(value_to_char2[1], 'C');
    EXPECT_EQ(value_to_char2[2], 'G');
    EXPECT_EQ(value_to_char2[3], 'T');
    EXPECT_EQ(value_to_char2[4], '-');

    constexpr std::array value_to_char3 = detail_union_composition<dna4, gap, dna5>::value_to_char_table();
    EXPECT_EQ(value_to_char3.size(), 10u);
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

    constexpr std::array value_to_char4 = detail_union_composition<dna5, gap, dna4>::value_to_char_table();
    EXPECT_EQ(value_to_char4.size(), 10u);
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

TEST(union_composition_detail_test, union_composition_char_to_value_table)
{
    constexpr std::array char_to_value1 = detail_union_composition<dna4>::char_to_value_table();
    EXPECT_EQ(char_to_value1.size(), 256u);
    EXPECT_EQ(char_to_value1['A'], 0);
    EXPECT_EQ(char_to_value1['C'], 1);
    EXPECT_EQ(char_to_value1['G'], 2);
    EXPECT_EQ(char_to_value1['T'], 3);

    constexpr std::array char_to_value2 = detail_union_composition<dna4, gap>::char_to_value_table();
    EXPECT_EQ(char_to_value2.size(), 256u);
    EXPECT_EQ(char_to_value2['A'], 0);
    EXPECT_EQ(char_to_value2['C'], 1);
    EXPECT_EQ(char_to_value2['G'], 2);
    EXPECT_EQ(char_to_value2['T'], 3);
    EXPECT_EQ(char_to_value2['-'], 4);

    constexpr std::array char_to_value3 = detail_union_composition<dna4, gap, dna5>::char_to_value_table();
    EXPECT_EQ(char_to_value3.size(), 256u);
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

    constexpr std::array char_to_value4 = detail_union_composition<dna5, gap, dna4>::char_to_value_table();
    EXPECT_EQ(char_to_value4.size(), 256u);
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
