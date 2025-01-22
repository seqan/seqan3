// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/alphabet/composite/alphabet_variant.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

template <typename... alphabet_types>
class detail_alphabet_variant : public seqan3::alphabet_variant<alphabet_types...>
{
public:
    using seqan3::alphabet_variant<alphabet_types...>::partial_sum_sizes;
    using seqan3::alphabet_variant<alphabet_types...>::rank_to_char_table;
    using seqan3::alphabet_variant<alphabet_types...>::char_to_rank_table;
    using seqan3::alphabet_variant<alphabet_types...>::first_valid_char_table;
};

TEST(alphabet_variant_detail_test, partial_sum_sizes)
{
    constexpr std::array partial_sum2 = detail_alphabet_variant<seqan3::dna4, seqan3::gap>::partial_sum_sizes;
    EXPECT_EQ(partial_sum2.size(), 3u);
    EXPECT_EQ(partial_sum2[0], 0);
    EXPECT_EQ(partial_sum2[1], 4);
    EXPECT_EQ(partial_sum2[2], 5);

    constexpr std::array partial_sum3 =
        detail_alphabet_variant<seqan3::dna4, seqan3::gap, seqan3::dna5>::partial_sum_sizes;
    EXPECT_EQ(partial_sum3.size(), 4u);
    EXPECT_EQ(partial_sum3[0], 0);
    EXPECT_EQ(partial_sum3[1], 4);
    EXPECT_EQ(partial_sum3[2], 5);
    EXPECT_EQ(partial_sum3[3], 10);

    constexpr std::array partial_sum4 =
        detail_alphabet_variant<seqan3::dna5, seqan3::gap, seqan3::dna4>::partial_sum_sizes;
    EXPECT_EQ(partial_sum4.size(), 4u);
    EXPECT_EQ(partial_sum4[0], 0);
    EXPECT_EQ(partial_sum4[1], 5);
    EXPECT_EQ(partial_sum4[2], 6);
    EXPECT_EQ(partial_sum4[3], 10);
}

TEST(alphabet_variant_detail_test, rank_to_char_table)
{
    constexpr std::array rank_to_char2 = detail_alphabet_variant<seqan3::dna4, seqan3::gap>::rank_to_char_table;
    EXPECT_EQ(rank_to_char2.size(), 5u);
    EXPECT_EQ(rank_to_char2[0], 'A');
    EXPECT_EQ(rank_to_char2[1], 'C');
    EXPECT_EQ(rank_to_char2[2], 'G');
    EXPECT_EQ(rank_to_char2[3], 'T');
    EXPECT_EQ(rank_to_char2[4], '-');

    constexpr std::array rank_to_char3 =
        detail_alphabet_variant<seqan3::dna4, seqan3::gap, seqan3::dna5>::rank_to_char_table;
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

    constexpr std::array rank_to_char4 =
        detail_alphabet_variant<seqan3::dna5, seqan3::gap, seqan3::dna4>::rank_to_char_table;
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

TEST(alphabet_variant_detail_test, char_to_rank_table)
{
    constexpr std::array char_to_rank2 = detail_alphabet_variant<seqan3::dna4, seqan3::gap>::char_to_rank_table;
    EXPECT_EQ(char_to_rank2.size(), 256u);
    EXPECT_EQ(char_to_rank2['A'], 0);
    EXPECT_EQ(char_to_rank2['C'], 1);
    EXPECT_EQ(char_to_rank2['G'], 2);
    EXPECT_EQ(char_to_rank2['T'], 3);
    EXPECT_EQ(char_to_rank2['-'], 4);

    constexpr std::array char_to_rank3 =
        detail_alphabet_variant<seqan3::dna4, seqan3::gap, seqan3::dna5>::char_to_rank_table;
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

    constexpr std::array char_to_rank4 =
        detail_alphabet_variant<seqan3::dna5, seqan3::gap, seqan3::dna4>::char_to_rank_table;
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

TEST(alphabet_variant_detail_test, first_valid_char_table)
{
    constexpr std::array first_valid_char_table2 =
        detail_alphabet_variant<seqan3::dna4, seqan3::gap>::first_valid_char_table;
    EXPECT_EQ(first_valid_char_table2.size(), 256u);
    EXPECT_EQ(first_valid_char_table2['A'], 0);
    EXPECT_EQ(first_valid_char_table2['C'], 0);
    EXPECT_EQ(first_valid_char_table2['G'], 0);
    EXPECT_EQ(first_valid_char_table2['T'], 0);
    EXPECT_EQ(first_valid_char_table2['-'], 1);
    EXPECT_EQ(first_valid_char_table2['X'], 2); // invalid char

    constexpr std::array first_valid_char_table3 =
        detail_alphabet_variant<seqan3::dna4, seqan3::gap, seqan3::dna5>::first_valid_char_table;
    EXPECT_EQ(first_valid_char_table3.size(), 256u);
    EXPECT_EQ(first_valid_char_table3['A'], 0);
    EXPECT_EQ(first_valid_char_table3['C'], 0);
    EXPECT_EQ(first_valid_char_table3['G'], 0);
    EXPECT_EQ(first_valid_char_table3['T'], 0);
    EXPECT_EQ(first_valid_char_table3['-'], 1);
    EXPECT_EQ(first_valid_char_table3['A'], 0);
    EXPECT_EQ(first_valid_char_table3['C'], 0);
    EXPECT_EQ(first_valid_char_table3['G'], 0);
    EXPECT_EQ(first_valid_char_table3['N'], 2);
    EXPECT_EQ(first_valid_char_table3['T'], 0);
    EXPECT_EQ(first_valid_char_table3['X'], 3); // invalid char

    constexpr std::array first_valid_char_table4 =
        detail_alphabet_variant<seqan3::dna5, seqan3::gap, seqan3::dna4>::first_valid_char_table;
    EXPECT_EQ(first_valid_char_table4.size(), 256u);
    EXPECT_EQ(first_valid_char_table4['A'], 0);
    EXPECT_EQ(first_valid_char_table4['C'], 0);
    EXPECT_EQ(first_valid_char_table4['G'], 0);
    EXPECT_EQ(first_valid_char_table4['N'], 0);
    EXPECT_EQ(first_valid_char_table4['T'], 0);
    EXPECT_EQ(first_valid_char_table4['-'], 1);
    EXPECT_EQ(first_valid_char_table4['A'], 0);
    EXPECT_EQ(first_valid_char_table4['C'], 0);
    EXPECT_EQ(first_valid_char_table4['G'], 0);
    EXPECT_EQ(first_valid_char_table4['T'], 0);
    EXPECT_EQ(first_valid_char_table4['X'], 3); // invalid char
}
