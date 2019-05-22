// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/view/convert.hpp>

#include "../alignment/aligned_sequence_test_template.hpp"

using namespace seqan3;

template <typename container_type>
    requires AlignedSequence<container_type>
class aligned_sequence<container_type> : public ::testing::Test
{
public:
    // Initializer function is needed for the typed test because the gapped_decorator
    // will be initialized differently than the naive vector<gapped<dna>>.
    void initialise_typed_test_container(container_type & container, dna4_vector const & target)
    {
        container.clear();
        for (auto & val : target)
        {
            container.push_back(assign_char_to(to_char(val), typename container_type::value_type{}));
        }
    }
};

using test_types = ::testing::Types<std::vector<gapped<dna4>>,
                                    std::vector<gapped<qualified<dna4, phred42>>>>;

INSTANTIATE_TYPED_TEST_CASE_P(container_of_gapped_alphabets, aligned_sequence, test_types);

TEST(aligned_sequence_debug_stream, multi_without_gaps)
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

    std::tuple<std::vector<gapped<dna4>>, std::vector<gapped<dna4>>, std::vector<gapped<dna4>>> const alignment
    {
        "GCGGGTCACTGAGGGCTGGGATGAGGACGGCCACCACTTCGAGGAGTCCCTTCACTACGAGGGCAGGGCCGTGGACATCACCACGTCAGACAGGGACAAGAGCAAGTA"
        "CGGCACCCTGTCCAGACTGGCGGTGGAAGCTG"_dna4 | view::persist | view::convert<gapped<dna4>>,
        "CTACGGCAGAAGAAGACATCCGAAAAAGCTGACACCTCTCGCCTACAAGCAGTTCATACCTAATGTCGCGGAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAA"
        "GATAACGCGCAATTCGGAGAGATTTAAAGAAC"_dna4 | view::persist | view::convert<gapped<dna4>>,
        "CTACGGCAGAAGAAGACATCCCAAGAAGCTGACACCTCTCGCCTACAAGCAGTTTATACCTAATGTCGCGGAGAAGACCTTAGGGGCCAGCGGCAGATACGAGGGCAA"
        "GATCACGCGCAATTCGGAGAGATTTAAAGAAC"_dna4 | view::persist | view::convert<gapped<dna4>>
    };

    std::ostringstream oss;
    debug_stream_type stream{oss};
    stream << alignment;
    EXPECT_EQ(expected, oss.str());
}

TEST(aligned_sequence_debug_stream, pair_with_gaps)
{
    std::string const expected
    {
        "      0     . \n"
        "        CUUC-G\n"
        "        ||   |\n"
        "        CU-NGG\n"
    };

    std::pair<std::vector<gapped<rna5>>, std::vector<gapped<rna5>>> const alignment
    {
        {'C'_rna5, 'U'_rna5, 'U'_rna5, 'C'_rna5, gap{}, 'G'_rna5},
        {'C'_rna5, 'U'_rna5, gap{}, 'N'_rna5, 'G'_rna5, 'G'_rna5}
    };

    std::ostringstream oss;
    debug_stream_type stream{oss};
    stream << alignment;
    EXPECT_EQ(expected, oss.str());
}
