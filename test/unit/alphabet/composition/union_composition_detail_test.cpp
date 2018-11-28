// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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
    using union_composition<alphabet_types...>::value_to_char;
    using union_composition<alphabet_types...>::char_to_value;
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

TEST(union_composition_detail_test, value_to_char)
{
    constexpr std::array value_to_char2 = detail_union_composition<dna4, gap>::value_to_char;
    EXPECT_EQ(value_to_char2.size(), 5u);
    EXPECT_EQ(value_to_char2[0], 'A');
    EXPECT_EQ(value_to_char2[1], 'C');
    EXPECT_EQ(value_to_char2[2], 'G');
    EXPECT_EQ(value_to_char2[3], 'T');
    EXPECT_EQ(value_to_char2[4], '-');

    constexpr std::array value_to_char3 = detail_union_composition<dna4, gap, dna5>::value_to_char;
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

    constexpr std::array value_to_char4 = detail_union_composition<dna5, gap, dna4>::value_to_char;
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

TEST(union_composition_detail_test, char_to_value)
{
    constexpr std::array char_to_value2 = detail_union_composition<dna4, gap>::char_to_value;
    EXPECT_EQ(char_to_value2.size(), 256u);
    EXPECT_EQ(char_to_value2['A'], 0);
    EXPECT_EQ(char_to_value2['C'], 1);
    EXPECT_EQ(char_to_value2['G'], 2);
    EXPECT_EQ(char_to_value2['T'], 3);
    EXPECT_EQ(char_to_value2['-'], 4);

    constexpr std::array char_to_value3 = detail_union_composition<dna4, gap, dna5>::char_to_value;
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

    constexpr std::array char_to_value4 = detail_union_composition<dna5, gap, dna4>::char_to_value;
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
