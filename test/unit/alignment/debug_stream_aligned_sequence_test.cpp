// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/range/views/convert.hpp>
#include <seqan3/range/views/to.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_rna5;

TEST(debug_stream_test, aligned_sequence_multi_without_gaps)
{
    std::string const expected
    {
        "      0     .    :    .    :    .    :    .    :    .    :\n"
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
        "        GAGGGCAAGATCACGCGCAATTCGGAGAGATTTAAAGAAC\n"
    };

    std::tuple<std::vector<seqan3::gapped<seqan3::dna4>>,
               std::vector<seqan3::gapped<seqan3::dna4>>,
               std::vector<seqan3::gapped<seqan3::dna4>>> const alignment
    {
        "GCGGGTCACTGAGGGCTGGGATGAGGACGGCCACCACTTCGAGGAGTCCCTTCACTACGAGGGCAGGGCCGTGGACATCACCACGTCAGACAGGGACAAGAGCAAGTA"
        "CGGCACCCTGTCCAGACTGGCGGTGGAAGCTG"_dna4 | seqan3::views::to<std::vector<seqan3::gapped<seqan3::dna4>>>,
        "CTACGGCAGAAGAAGACATCCGAAAAAGCTGACACCTCTCGCCTACAAGCAGTTCATACCTAATGTCGCGGAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAA"
        "GATAACGCGCAATTCGGAGAGATTTAAAGAAC"_dna4 | seqan3::views::to<std::vector<seqan3::gapped<seqan3::dna4>>>,
        "CTACGGCAGAAGAAGACATCCCAAGAAGCTGACACCTCTCGCCTACAAGCAGTTTATACCTAATGTCGCGGAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAA"
        "GATCACGCGCAATTCGGAGAGATTTAAAGAAC"_dna4 | seqan3::views::to<std::vector<seqan3::gapped<seqan3::dna4>>>
    };

    std::ostringstream oss;
    seqan3::debug_stream_type stream{oss};
    stream << alignment;
    EXPECT_EQ(expected, oss.str());
}

TEST(debug_stream_test, aligned_sequence_pair_with_gaps)
{
    std::string const expected
    {
        "      0     . \n"
        "        CUUC-G\n"
        "        ||   |\n"
        "        CU-NGG\n"
    };

    std::pair<std::vector<seqan3::gapped<seqan3::rna5>>, std::vector<seqan3::gapped<seqan3::rna5>>> const alignment
    {
        {'C'_rna5, 'U'_rna5, 'U'_rna5, 'C'_rna5, seqan3::gap{}, 'G'_rna5},
        {'C'_rna5, 'U'_rna5, seqan3::gap{}, 'N'_rna5, 'G'_rna5, 'G'_rna5}
    };

    std::ostringstream oss;
    seqan3::debug_stream_type stream{oss};
    stream << alignment;
    EXPECT_EQ(expected, oss.str());
}
