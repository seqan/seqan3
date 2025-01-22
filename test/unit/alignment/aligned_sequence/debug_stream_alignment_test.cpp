// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/aligned_sequence/debug_stream_alignment.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/utility/range/to.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_rna5;

TEST(debug_stream_test, multiple_alignment_without_gaps)
{
    std::string const expected{"      0     .    :    .    :    .    :    .    :    .    :\n"
                               "        GCGGGTCACTGAGGGCTGGGATGAGGACGGCCACCACTTCGAGGAGTCCC\n"
                               "            | ||      |        |  |       |   |||   |    |\n"
                               "        CTACGGCAGAAGAAGACATCCGAAAAAGCTGACACCTCTCGCCTACAAGC\n"
                               "        ||||||||||||||||||||| || |||||||||||||||||||||||||\n"
                               "        CTACGGCAGAAGAAGACATCCCAAGAAGCTGACACCTCTCGCCTACAAGC\n"
                               "\n"
                               "     50     .    :    .    :    .    :    .    :    .    :\n"
                               "        TTCACTACGAGGGCAGGGCCGTGGACATCACCACGTCAGACAGGGACAAG\n"
                               "            |            || | | | | |     | |   | |     | \n"
                               "        AGTTCATACCTAATGTCGCGGAGAAGACCTTAGGGGCCAGCGGCAGATAC\n"
                               "        |||| |||||||||||||||||||||||||||||||||||||||||||||\n"
                               "        AGTTTATACCTAATGTCGCGGAGAAGACCTTAGGGGCCAGCGGCAGATAC\n"
                               "\n"
                               "    100     .    :    .    :    .    :    .    :\n"
                               "        AGCAAGTACGGCACCCTGTCCAGACTGGCGGTGGAAGCTG\n"
                               "               |    || |          |    |  |||   \n"
                               "        GAGGGCAAGATAACGCGCAATTCGGAGAGATTTAAAGAAC\n"
                               "        ||||||||||| ||||||||||||||||||||||||||||\n"
                               "        GAGGGCAAGATCACGCGCAATTCGGAGAGATTTAAAGAAC\n"};

    std::tuple<std::vector<seqan3::gapped<seqan3::dna4>>,
               std::vector<seqan3::gapped<seqan3::dna4>>,
               std::vector<seqan3::gapped<seqan3::dna4>>> const alignment{
        "GCGGGTCACTGAGGGCTGGGATGAGGACGGCCACCACTTCGAGGAGTCCCTTCACTACGAGGGCAGGGCCGTGGACATCACCACGTCAGACAGGGACAAGAGCAAGTA"
        "CGGCACCCTGTCCAGACTGGCGGTGGAAGCTG"_dna4
            | seqan3::ranges::to<std::vector<seqan3::gapped<seqan3::dna4>>>(),
        "CTACGGCAGAAGAAGACATCCGAAAAAGCTGACACCTCTCGCCTACAAGCAGTTCATACCTAATGTCGCGGAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAA"
        "GATAACGCGCAATTCGGAGAGATTTAAAGAAC"_dna4
            | seqan3::ranges::to<std::vector<seqan3::gapped<seqan3::dna4>>>(),
        "CTACGGCAGAAGAAGACATCCCAAGAAGCTGACACCTCTCGCCTACAAGCAGTTTATACCTAATGTCGCGGAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAA"
        "GATCACGCGCAATTCGGAGAGATTTAAAGAAC"_dna4
            | seqan3::ranges::to<std::vector<seqan3::gapped<seqan3::dna4>>>()};

    std::ostringstream oss;
    seqan3::debug_stream_type stream{oss};
    stream << alignment;
    EXPECT_EQ(expected, oss.str());
}

TEST(debug_stream_test, pairwise_alignment_with_gaps)
{
    std::string const expected{"      0     . \n"
                               "        CUUC-G\n"
                               "        ||   |\n"
                               "        CU-NGG\n"};

    std::pair<std::vector<seqan3::gapped<seqan3::rna5>>, std::vector<seqan3::gapped<seqan3::rna5>>> const alignment{
        {'C'_rna5, 'U'_rna5, 'U'_rna5, 'C'_rna5, seqan3::gap{}, 'G'_rna5},
        {'C'_rna5, 'U'_rna5, seqan3::gap{}, 'N'_rna5, 'G'_rna5, 'G'_rna5}};

    std::ostringstream oss;
    seqan3::debug_stream_type stream{oss};
    stream << alignment;
    EXPECT_EQ(expected, oss.str());
}
