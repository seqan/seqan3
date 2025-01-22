// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>
#include <sstream>
#include <string>

#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sequence_file/format_fastq.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/utility/views/convert.hpp>

#include "sequence_file_format_test_template.hpp"

using seqan3::operator""_dna5;
using seqan3::operator""_phred42;

template <>
struct sequence_file_read<seqan3::format_fastq> : public sequence_file_data
{
    std::string standard_input{"@ID1\n"
                               "ACGTTTTTTTTTTTTTTT\n"
                               "+\n"
                               "!##$%&'()*+,-./++-\n"
                               "@ID2\n"
                               "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
                               "+\n"
                               "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
                               "@ID3 lala\n"
                               "ACGTTTA\n"
                               "+\n"
                               "!!!!!!!\n"};

    std::string illegal_alphabet_character_input{"@ID1\n"
                                                 "ACGTTPTTTTTTTTTTTT\n"
                                                 "+\n"
                                                 "!##$%&'()*+,-./++-\n"};

    std::string standard_output{"@ID1\n"
                                "ACGTTTTTTTTTTTTTTT\n"
                                "+\n"
                                "!##$%&'()*+,-./++-\n"
                                "@ID2\n"
                                "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
                                "+\n"
                                "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
                                "@ID3 lala\n"
                                "ACGTTTA\n"
                                "+\n"
                                "!!!!!!!\n"};

    std::string no_or_ill_formatted_id_input{"#ID1\n"
                                             "ACGTTPTTTTTTTTTTTT\n"
                                             "+\n"
                                             "!##$%&'()*+,-./++-\n"};
};

// ---------------------------------------------------------------------------------------------------------------------
// parametrized tests
// ---------------------------------------------------------------------------------------------------------------------

INSTANTIATE_TYPED_TEST_SUITE_P(fastq, sequence_file_read, seqan3::format_fastq, );
INSTANTIATE_TYPED_TEST_SUITE_P(fastq, sequence_file_write, seqan3::format_fastq, );

// ----------------------------------------------------------------------------
// reading
// ----------------------------------------------------------------------------

struct read : public sequence_file_data
{
    std::string input{"@ID1\n"
                      "ACGTTTTTTTTTTTTTTT\n"
                      "+\n"
                      "!##$%&'()*+,-./++-\n"
                      "@ID2\n"
                      "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
                      "+\n"
                      "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
                      "@ID3 lala\n"
                      "ACGTTTA\n"
                      "+\n"
                      "!!!!!!!\n"};

    seqan3::sequence_file_input_options<seqan3::dna15> options{};

    void do_read_test(std::string const & input)
    {
        std::stringstream istream{input};
        seqan3::sequence_file_input fin{istream, seqan3::format_fastq{}};
        fin.options = options;

        auto it = fin.begin();
        for (unsigned i = 0; i < 3; ++i, ++it)
        {
            EXPECT_RANGE_EQ((*it).id(), ids[i]);
            EXPECT_RANGE_EQ((*it).sequence(), seqs[i]);
            EXPECT_RANGE_EQ((*it).base_qualities(), quals[i]);
        }
    }
};

TEST_F(read, newline_before_eof)
{
    input = "@ID1\n"
            "ACGTTTTTTTTTTTTTTT\n"
            "+\n"
            "!##$%&'()*+,-./++-\n"
            "@ID2\n"
            "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
            "+\n"
            "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
            "@ID3 lala\n"
            "ACGTTTA\n"
            "+\n"
            "!!!!!!!";
    do_read_test(input);
}

TEST_F(read, whitespace_in_seq_qual)
{
    input = "@ID1\n"
            "ACGTTTTTTTT\nTTTTTTT\n"
            "+\n"
            "!##$\n%&'()*+,-./++-\n"
            "@ID2\n"
            "ACGTTTTTTTTTT\r\nTTTTTT\r\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\r\nTTTTTTTTTTTTTTTTTTTTTTT\n"
            "+\n"
            "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBD\nEBDEB\nDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
            "@ID3 lala\n"
            "ACGTT\nTA\n"
            "+\n"
            "!!!!!\n!!\n";
    do_read_test(input);
}

TEST_F(read, double_id_style)
{
    input = "@ID1\n"
            "ACGTTTTTTTTTTTTTTT\n"
            "+ID1\n"
            "!##$%&'()*+,-./++-\n"
            "@ID2\n"
            "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
            "+ID2\n"
            "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
            "@ID3 lala\n"
            "ACGTTTA\n"
            "+ID3 lala\n"
            "!!!!!!!\n";
    do_read_test(input);
}

TEST_F(read, mixed_issues)
{
    input = "@ID1\n"
            "ACGTTTTTTTT\nTTTTTTT\n"
            "+\n"
            "!##$\n%&'()*+,-./++-\n"
            "@ID2\n"
            "ACGTTTTTTTTTT\r\nTTTTTT\r\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\r\nTTTTTTTTTTTTTTTTTTTTTTT\n"
            "+\n"
            "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBD\nEBDEB\nDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
            "@ID3 lala\n"
            "ACGTT\nTA\n"
            "+ID3 lala\n"
            "!!!!!\n!!";
    do_read_test(input);
}

TEST_F(read, only_qual)
{
    std::stringstream istream{input};
    seqan3::sequence_file_input fin{istream, seqan3::format_fastq{}, seqan3::fields<seqan3::field::qual>{}};

    auto it = fin.begin();
    for (unsigned i = 0; i < 3; ++i, ++it)
        EXPECT_RANGE_EQ(std::get<0>(*it), quals[i]);
}

//TODO fail_no_2nd_id
TEST_F(read, fail_no_seq_after_id)
{
    std::stringstream istream{input = "@ID1\n"
                                      "ACGTTTTTTTTTTTTTTT"};

    seqan3::sequence_file_input fin{istream, seqan3::format_fastq{}};
    EXPECT_THROW(fin.begin(), seqan3::unexpected_end_of_input);
}

//TODO fail_no_quals
//TODO fail_quals_shorter_seq

// ----------------------------------------------------------------------------
// writing
// ----------------------------------------------------------------------------

struct write : public ::testing::Test
{
    std::vector<seqan3::dna5_vector> seqs{
        "ACGT"_dna5,
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN"_dna5,
        "GGAGTATAATATATATATATAT"_dna5};

    std::vector<std::string> ids{"TEST 1", "Test2", "Test3"};

    std::vector<std::vector<seqan3::phred42>> quals{
        {"!##$"_phred42},
        {"!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE"_phred42},
        {"!!*+,-./+*+,-./+!!FF!!"_phred42},
    };

    seqan3::sequence_file_output_options options{};

    std::ostringstream ostream;

    void do_write_test()
    {
        seqan3::sequence_file_output fout{ostream, seqan3::format_fastq{}};
        fout.options = options;

        for (unsigned i = 0; i < 3; ++i)
            EXPECT_NO_THROW((fout.emplace_back(seqs[i], ids[i], quals[i])));

        ostream.flush();
    }
};

TEST_F(write, arg_handling_qual_missing)
{
    seqan3::sequence_file_output fout{ostream,
                                      seqan3::format_fastq{},
                                      seqan3::fields<seqan3::field::id, seqan3::field::seq>{}};
    EXPECT_THROW((fout.emplace_back(ids[0], seqs[0])), std::logic_error);
}

TEST_F(write, arg_handling_qual_empty)
{
    seqan3::sequence_file_output fout{ostream, seqan3::format_fastq{}};
    EXPECT_THROW((fout.emplace_back(seqs[0], ids[0], std::string_view{""})), std::runtime_error);
}

TEST_F(write, options_fastq_double_id)
{
    options.fastq_double_id = true;

    std::string comp{"@TEST 1\n"
                     "ACGT\n"
                     "+TEST 1\n"
                     "!##$\n"
                     "@Test2\n"
                     "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN\n"
                     "+Test2\n"
                     "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE\n"
                     "@Test3\n"
                     "GGAGTATAATATATATATATAT\n"
                     "+Test3\n"
                     "!!*+,-./+*+,-./+!!FF!!\n"};

    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, options_add_carriage_return)
{
    options.add_carriage_return = true;

    std::string comp{"@TEST 1\r\n"
                     "ACGT\r\n"
                     "+\r\n"
                     "!##$\r\n"
                     "@Test2\r\n"
                     "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN\r\n"
                     "+\r\n"
                     "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE\r\n"
                     "@Test3\r\n"
                     "GGAGTATAATATATATATATAT\r\n"
                     "+\r\n"
                     "!!*+,-./+*+,-./+!!FF!!\r\n"};

    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, options_all)
{
    options.add_carriage_return = true;
    options.fastq_double_id = true;

    std::string comp{"@TEST 1\r\n"
                     "ACGT\r\n"
                     "+TEST 1\r\n"
                     "!##$\r\n"
                     "@Test2\r\n"
                     "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN\r\n"
                     "+Test2\r\n"
                     "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE\r\n"
                     "@Test3\r\n"
                     "GGAGTATAATATATATATATAT\r\n"
                     "+Test3\r\n"
                     "!!*+,-./+*+,-./+!!FF!!\r\n"};
    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}
