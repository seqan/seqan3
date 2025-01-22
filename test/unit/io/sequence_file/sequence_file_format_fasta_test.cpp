// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>
#include <sstream>
#include <string>

#include <seqan3/alphabet/nucleotide/dna15.hpp>

#include "sequence_file_format_test_template.hpp"

template <>
struct sequence_file_read<seqan3::format_fasta> : public sequence_file_data
{
    std::string standard_input{">ID1\n"
                               "ACGTTTTTTTTTTTTTTT\n"
                               ">ID2\n"
                               "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
                               ">ID3 lala\n"
                               "ACGTTTA\n"};

    std::string illegal_alphabet_character_input{
        ">ID1\n"
        "ACGPTTTTTTTTTTTTTT\n"
        ">ID2\n"
        "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        ">ID3 lala\n"
        "ACGTTTA\n"};

    std::string standard_output{">ID1\n"
                                "ACGTTTTTTTTTTTTTTT\n"
                                ">ID2\n"
                                "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\nTT\n"
                                //                                             linebreak inserted after 80 char  ^
                                ">ID3 lala\n"
                                "ACGTTTA\n"};

    std::string no_or_ill_formatted_id_input{"! ID1\n"
                                             "ACGTTTTTTTTTTTTTTT\n"};
};

// ---------------------------------------------------------------------------------------------------------------------
// parametrized tests
// ---------------------------------------------------------------------------------------------------------------------

INSTANTIATE_TYPED_TEST_SUITE_P(fasta, sequence_file_read, seqan3::format_fasta, );
INSTANTIATE_TYPED_TEST_SUITE_P(fasta, sequence_file_write, seqan3::format_fasta, );

// ----------------------------------------------------------------------------
// reading
// ----------------------------------------------------------------------------

struct read : public sequence_file_data
{
    seqan3::sequence_file_input_options<seqan3::dna15> options{};

    std::string id;
    seqan3::dna5_vector seq;

    void do_read_test(std::string const & input)
    {
        std::stringstream istream{input};

        seqan3::sequence_file_input fin{istream,
                                        seqan3::format_fasta{},
                                        seqan3::fields<seqan3::field::id, seqan3::field::seq>{}};
        fin.options = options;

        auto it = fin.begin();
        for (unsigned i = 0; i < 3; ++i, ++it)
        {
            EXPECT_EQ((*it).id(), ids[i]);
            EXPECT_RANGE_EQ((*it).sequence(), seqs[i]);
        }
    }
};

TEST_F(read, newline_before_eof)
{
    std::string input{">ID1\n"
                      "ACGTTTTTTTTTTTTTTT\n"
                      ">ID2\n"
                      "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
                      ">ID3 lala\n"
                      "ACGTTTA"};
    do_read_test(input);
}

TEST_F(read, noblank_before_id)
{
    std::string input{">ID1\n"
                      "ACGTTTTTTTTTTTTTTT\n"
                      ">ID2\n"
                      "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
                      ">ID3 lala\n"
                      "ACGTTTA\n"};
    do_read_test(input);
}

TEST_F(read, whitespace_in_seq)
{
    std::string input{">ID1\n"
                      "ACGTTTT\n\nTTTTTTTTTTT\n"
                      "\n"
                      ">ID2\n"
                      "ACGTTTT\t\tTTTTTTTTTTT\t\nTTTTTTTTTTT\vTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\rTTTTTTTTTTTTTTTTT\n"
                      ">ID3 lala\n"
                      "ACGT\fTTA\n"};
    do_read_test(input);
}

struct char_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = char;
    using sequence_legal_alphabet = seqan3::dna4;
};
using sequence_file_type = seqan3::sequence_file_input<char_traits,
                                                       seqan3::fields<seqan3::field::id, seqan3::field::seq>,
                                                       seqan3::type_list<seqan3::format_fasta>>;

TEST_F(read, whitespace_in_seq_char_alphabet)
{
    std::string input{">ID1\n"
                      "ACGTTTT\n\nTTTTTTTTTTT\n"
                      "\n"
                      ">ID2\n"
                      "ACGTTTT\t\tTTTTTTTTTTT\t\nTTTTTTTTTTT\vTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\rTTTTTTTTTTTTTTTTT\n"
                      ">ID3 lala\n"
                      "ACGT\fTTA\n"};
    std::stringstream istream{input};

    sequence_file_type fin{istream, seqan3::format_fasta{}};

    auto it = fin.begin();
    for (unsigned i = 0; i < 3; ++i, ++it)
    {
        EXPECT_EQ((*it).id(), ids[i]);
        EXPECT_RANGE_EQ((*it).sequence(), seqs[i] | seqan3::views::to_char);
    }
}

TEST_F(read, digits_in_seq)
{
    std::string input{">ID1\n"
                      "10  ACGTTTTTTTTTTTTTTT\n"
                      ">ID2\n"
                      "  80 ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT  900"
                      "1000 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
                      ">ID3 lala\n"
                      "ACGT9T5T2A\n"};
    do_read_test(input);
}

TEST_F(read, old_id_style)
{
    std::string input{"; ID1\n"
                      "ACGTTTTTTTTTTTTTTT\n"
                      "; ID2\n"
                      "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
                      "; ID3 lala\n"
                      "ACGTTTA\n"};

    do_read_test(input);
}

TEST_F(read, mixed_issues)
{
    std::string input{
        ">ID1\n"
        "ACGTTTT\n\nTTTTTTTTTTT\n"
        "\n"
        ";ID2\n"
        "ACGTTTT\t75\tTTTTTTTTTTT\t\nTTTTTTTTTTT9\vTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\rTTTTTTTTTTTTTTTTT\n"
        ">ID3 lala\n"
        "ACGT\fTTA"};
    do_read_test(input);
}

TEST_F(read, fail_no_newline_after_id)
{
    std::string input{">ID1"
                      "ACGTTTTTTTTTTTTTTT"};
    std::stringstream istream{input};

    seqan3::sequence_file_input fin{istream, seqan3::format_fasta{}};
    EXPECT_THROW(fin.begin(), seqan3::unexpected_end_of_input);
}

TEST_F(read, fail_no_newline_after_truncate_id)
{
    std::string input{">ID1 to_be_truncated"
                      "ACGTTTTTTTTTTTTTTT"};
    std::stringstream istream{input};

    seqan3::sequence_file_input fin{istream, seqan3::format_fasta{}};
    fin.options.truncate_ids = true;
    EXPECT_THROW(fin.begin(), seqan3::unexpected_end_of_input);
}

// ----------------------------------------------------------------------------
// writing
// ----------------------------------------------------------------------------

struct write : public ::testing::Test
{
    std::vector<seqan3::dna5_vector> seqs{
        "ACGT"_dna5,
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN"_dna5,
        "GGAGTATAATATATATATATATAT"_dna5};

    std::vector<std::string> ids{"TEST 1", "Test2", "Test3"};

    seqan3::sequence_file_output_options options{};

    std::ostringstream ostream;

    void do_write_test()
    {
        seqan3::sequence_file_output fout{ostream,
                                          seqan3::format_fasta{},
                                          seqan3::fields<seqan3::field::seq, seqan3::field::id>{}};
        fout.options = options;

        for (unsigned i = 0; i < 3; ++i)
            EXPECT_NO_THROW((fout.emplace_back(seqs[i], ids[i])));

        ostream.flush();
    }
};

TEST_F(write, options_letters_per_line)
{
    options.fasta_letters_per_line = 7;

    std::string comp{
        ">TEST 1\n"
        "ACGT\n"
        ">Test2\n"
        "AGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\n"
        "AGGCTGN\n"
        ">Test3\n"
        "GGAGTAT\nAATATAT\nATATATA\nTAT\n"};

    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, options_legacy_id_marker)
{
    options.fasta_legacy_id_marker = true;

    std::string comp{";TEST 1\n"
                     "ACGT\n"
                     ";Test2\n"
                     "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGG\nCTGNAGGCTGN\n"
                     //                                             linebreak inserted after 80 char  ^
                     ";Test3\n"
                     "GGAGTATAATATATATATATATAT\n"};

    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, options_blank_before_id)
{
    options.fasta_blank_before_id = true;

    std::string comp{"> TEST 1\n"
                     "ACGT\n"
                     "> Test2\n"
                     "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGG\nCTGNAGGCTGN\n"
                     //                                             linebreak inserted after 80 char  ^
                     "> Test3\n"
                     "GGAGTATAATATATATATATATAT\n"};

    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, options_add_carriage_return)
{
    options.add_carriage_return = true;

    std::string comp{
        ">TEST 1\r\n"
        "ACGT\r\n"
        ">Test2\r\n"
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGG\r\nCTGNAGGCTGN\r\n"
        //                                             linebreak inserted after 80 char  ^
        ">Test3\r\n"
        "GGAGTATAATATATATATATATAT\r\n"};

    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, options_all)
{
    options.add_carriage_return = true;
    options.fasta_blank_before_id = true;
    options.fasta_legacy_id_marker = true;
    options.fasta_letters_per_line = 21;

    std::string comp{
        "; TEST 1\r\n"
        "ACGT\r\n"
        "; Test2\r\n"
        "AGGCTGNAGGCTGNAGGCTGN\r\nAGGCTGNAGGCTGNAGGCTGN\r\nAGGCTGNAGGCTGNAGGCTGN\r\nAGGCTGNAGGCTGNAGGCTGN\r\n"
        "AGGCTGN\r\n"
        "; Test3\r\n"
        "GGAGTATAATATATATATATA\r\nTAT\r\n"};
    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}
