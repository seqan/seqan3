// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================

/*!\file
 * \brief Provides tests for the seqan3::alignment class.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#include <gtest/gtest.h>

#include <sstream>

#include <seqan3/align/alignment.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>

using namespace seqan3;
using namespace seqan3::literal;

TEST(alignment_class_test, constructor_and_ostream)
{
    alignment align("GCGGGTCACTGAGGGCTGGGATGAGGACGGCCACCACTTCGAGGAGTCCCTTCACTACGAGGGCAGGGCC"
                    "GTGGACATCACCACGTCAGACAGGGACAAGAGCAAGTACGGCACCCTGTCCAGACTGGCGGTGGAAGCTG"_dna4,
                    "CTACGGCAGAAGAAGACATCCGAAAAAGCTGACACCTCTCGCCTACAAGCAGTTCATACCTAATGTCGCG"
                    "GAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAAGATAACGCGCAATTCGGAGAGATTTAAAGAAC"_dna4,
                    "CTACGGCAGAAGAAGACATCCCAAGAAGCTGACACCTCTCGCCTACAAGCAGTTTATACCTAATGTCGCG"
                    "GAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAAGATCACGCGCAATTCGGAGAGATTTAAAGAAC"_dna4);
    std::stringstream stream;
    stream << align;

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

    EXPECT_EQ(expected, stream.str());
}

TEST(alignment_class_test, column_view)
{
    alignment align("GCGG"_dna4, "CTAC"_dna4, "CTAC"_dna4);
    column_view_type<dna4_vector, dna4_vector, dna4_vector> columns = column_view(align);
    auto iter = columns.begin();
    // 1st alignment iter
    EXPECT_EQ(std::get<0>(*iter), dna4{dna4::G});
    EXPECT_EQ(std::get<1>(*iter), dna4{dna4::C});
    EXPECT_EQ(std::get<2>(*iter), dna4{dna4::C});
    ++iter; // 2nd iter
    EXPECT_TRUE(iter > columns.begin());
    EXPECT_EQ(std::get<1>(*iter), dna4{dna4::T});
    iter += 2; // 4th iter
    EXPECT_EQ(std::get<1>(*iter), dna4{dna4::C});
    --iter; // 3rd iter
    EXPECT_EQ(std::get<1>(*iter), dna4{dna4::A});
    ++iter;
    iter++; // end
    EXPECT_TRUE(iter == columns.end());

    std::stringstream stream;
    std::for_each(columns.begin(), columns.end(), [&stream] (auto const & col)
    {
        stream << std::get<0>(col) << std::get<1>(col) << std::get<2>(col) << ' ';
    });
    EXPECT_EQ("GCC CTT GAA GCC ", stream.str());
}

TEST(alignment_class_test, column_view_deduced)
{
    alignment align("AUUGN"_rna5, "AUUGN"_rna5);
    EXPECT_EQ(align.depth, 2);
    for (auto const & col : column_view(align))
        EXPECT_EQ(std::get<0>(col), std::get<1>(col));
}

TEST(alignment_class_test, depth)
{
    alignment align("GCGG"_dna4, "CTAC"_dna4, "CTAC"_dna4);
    EXPECT_EQ(align.depth, 3);
}

TEST(alignment_class_test, different_sequence_types)
{
    std::string const expected
    {
        "      0     \n"
        "        CTTC\n"
        "        ||  \n"
        "        CTAN\n"
        "        | | \n"
        "        CUAC\n"
    };

    alignment align("CTTC"_dna4, "CTAN"_dna5, "CUAC"_rna4);
    EXPECT_EQ(align.depth, 3);

    std::stringstream stream;
    stream << align;
    EXPECT_EQ(stream.str(), expected);
}
