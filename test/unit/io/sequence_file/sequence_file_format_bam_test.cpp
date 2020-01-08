// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/alignment_file/input.hpp>
#include <seqan3/io/alignment_file/output.hpp>
#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/test/pretty_printing.hpp>

#include "sequence_file_format_test_template.hpp"

using namespace seqan3;

template <>
struct sequence_file_read<format_bam> : public sequence_file_data<format_bam>
{
    std::string standard_input
    {
        // BGZF HEADER
        '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43',
        '\x02', '\x00', '\x92', '\x00',
        // CDATA
        '\x73', '\x72', '\xF4', '\x65', '\xE4', '\x66', '\x60', '\x60', '\x70', '\xF0', '\x70', '\xE1', '\x0C', '\xF3',
        '\xB3', '\x32', '\xD4', '\x33', '\xE3', '\x02', '\xF2', '\x18', '\xEC', '\x81', '\xF8', '\x3F', '\x14', '\xB0',
        '\x30', '\x78', '\x08', '\x81', '\xC4', '\x40', '\x04', '\x4C', '\x0C', '\xC4', '\xF7', '\x74', '\x31', '\x64',
        '\x10', '\xF2', '\xE8', '\x80', '\x00', '\x06', '\x26', '\x26', '\x66', '\x16', '\x56', '\x36', '\x76', '\x0E',
        '\x4E', '\x2E', '\x6E', '\x1E', '\x5E', '\x3E', '\x2E', '\x2E', '\x9E', '\xF9', '\x48', '\xCA', '\x61', '\x46',
        '\x04', '\x21', '\x89', '\x41', '\x8C', '\x30', '\x42', '\x18', '\x01', '\x35', '\x09', '\x27', '\x05', '\xB2',
        '\x02', '\xC9', '\x06', '\xB0', '\x4D', '\xA8', '\x84', '\xA2', '\xB2', '\x0A', '\x49', '\xC8', '\x04', '\xC9',
        '\x39', '\x9C', '\x50', '\x5F', '\xB2', '\x23', '\x89', '\x41', '\x9C', '\x68', '\xAC', '\x90', '\x93', '\x98',
        '\x93', '\x08', '\x72', '\xA7', '\x00', '\x48', '\x00', '\x08', '\x00',
        // CRC32
        '\x90', '\xD9', '\x82', '\xA8',
        // ISIZE
        '\x35', '\x01', '\x00', '\x00',
        // EOF-marker
        '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43',
        '\x02', '\x00', '\x1B', '\x00', '\x03', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00'
    };

    std::string illegal_alphabet_character_input
    {
        // HEADER
        '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43',
        '\x02', '\x00', '\x21', '\x00',
        // CDATA
        '\x73', '\x72', '\xF4', '\x65', '\x64', '\x80', '\x02', '\x00', '\x7C', '\xB1', '\x74', '\xCC', '\x0C', '\x00',
        '\x00', '\x00', '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00',
        '\x42', '\x43', '\x02', '\x00', '\x4D', '\x00', '\xB3', '\x67', '\x60', '\x60', '\xF8', '\x0F', '\x05', '\x2C',
        '\x0C', '\x1E', '\x42', '\x0C', '\x0C', '\x2C', '\x0C', '\x42', '\x48', '\x62', '\x40', '\x26', '\x83', '\xA7',
        '\x8B', '\x21', '\x83', '\x90', '\x47', '\x47', '\x47', '\x87', '\x05', '\x10', '\x33', '\x30', '\x31', '\x31',
        '\xB3', '\xB0', '\xB2', '\xB1', '\x73', '\x70', '\x72', '\x71', '\xF3', '\xF0', '\xF2', '\x71', '\x71', '\xF1',
        '\x00', '\x00',
        // CRC32
        '\xCA', '\xAB', '\xDF', '\x63',
        // ISIZE
        '\x43', '\x00', '\x00', '\x00',
        // EOF-marker
        '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43',
        '\x02', '\x00', '\x1B', '\x00', '\x03', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00'
    };

    std::string standard_output {standard_input};

    std::string no_or_ill_formatted_id_input
    {
        '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43',
        '\x02', '\x00', '\x21', '\x00', '\x73', '\x72', '\xF4', '\x65', '\x64', '\x80', '\x02', '\x00', '\x7C', '\xB1',
        '\x74', '\xCC', '\x0C', '\x00', '\x00', '\x00', '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00',
        '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43', '\x02', '\x00', '\x48', '\x00', '\xB3', '\x65', '\x60', '\x60',
        '\xF8', '\x0F', '\x05', '\x4C', '\x0C', '\x1E', '\x42', '\x0C', '\x0C', '\x2C', '\x0C', '\x42', '\x48', '\x62',
        '\x40', '\x26', '\x83', '\x16', '\x83', '\x90', '\x47', '\x07', '\x04', '\x30', '\x30', '\x31', '\x31', '\xB3',
        '\xB0', '\xB2', '\xB1', '\x73', '\x70', '\x72', '\x71', '\xF3', '\xF0', '\xF2', '\x71', '\x71', '\xF1', '\x00',
        '\x00', '\x6B', '\x33', '\x7C', '\xAB', '\x41', '\x00', '\x00', '\x00', '\x1F', '\x8B', '\x08', '\x04', '\x00',
        '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43', '\x02', '\x00', '\x1B', '\x00', '\x03',
        '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00'
	};
};

// ---------------------------------------------------------------------------------------------------------------------
// parametrized tests
// ---------------------------------------------------------------------------------------------------------------------

INSTANTIATE_TYPED_TEST_CASE_P(bam, sequence_file_read, format_bam);
INSTANTIATE_TYPED_TEST_CASE_P(bam, sequence_file_write, format_bam);

// ----------------------------------------------------------------------------
// reading
// ----------------------------------------------------------------------------
struct read_bam : sequence_file_read<format_bam>
{
    sequence_file_input_options<dna5, false> options{};

    std::string single_record_bam_raw
    {
        // First comment is byte offset of the first byte.
        // Header - 23 byte
        /*00*/ '\x42', '\x41', '\x4D', '\x01',                                             // BAM1
        /*04*/ '\x0B', '\x00', '\x00', '\x00',                                             // header size (11)
        /*08*/ '\x40', '\x48', '\x44', '\x09', '\x56', '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', // @HD VN:1.6
        /*19*/ '\x00', '\x00', '\x00', '\x00',                                             // n_ref (0)
        // Alignment record - 67 byte
        /*23*/ '\x3F', '\x00', '\x00', '\x00',                                             // block_size (63)
        /*27*/ '\xFF', '\xFF', '\xFF', '\xFF',                                             // refID (-1)
        /*31*/ '\xFF', '\xFF', '\xFF', '\xFF',                                             // pos (-1)
        /*35*/ '\x04',                                                                     // l_read_name (4)
        /*36*/ '\x00',                                                                     // mapq (0)
        /*37*/ '\x48', '\x12',                                                             // bin (4680)
        /*39*/ '\x00', '\x00',                                                             // n_cigar_bin (0)
        /*41*/ '\x00', '\x00',                                                             // flag (0)
        /*43*/ '\x12', '\x00', '\x00', '\x00',                                             // l_seq (18)
        /*47*/ '\xFF', '\xFF', '\xFF', '\xFF',                                             // next_refId (-1)
        /*51*/ '\xFF', '\xFF', '\xFF', '\xFF',                                             // next_pos (-1)
        /*55*/ '\x00', '\x00', '\x00', '\x00',                                             // tlen (0)
        /*59*/ '\x49', '\x44', '\x31', '\x00',                                             // read_name (ID1)
                                                                                           // cigar (*)
        /*63*/ '\x12', '\x48', '\x88', '\x88', '\x88', '\x88', '\x88', '\x88', '\x88',     // seq  (ACGTTTTTTTTTTTTTTT)
        /*72*/ '\x00', '\x02', '\x02', '\x03', '\x04', '\x05', '\x06', '\x07', '\x08',     // qual (!##$%&'()*+,-./++-)
        /*81*/ '\x09', '\x0A', '\x0B', '\x0C', '\x0D', '\x0E', '\x0A', '\x0A', '\x0C'
    };

    void do_read_test(std::string const & input)
    {
        std::stringstream istream{input};
        std::stringstream buffer_stream{};

        // First convert the sam into a bam format.
        {
            alignment_file_input fin{istream, format_sam{}};
            alignment_file_output fout{buffer_stream, format_bam{}};
            fout = fin;
        }

        // Then read in bam as sequence file.
        sequence_file_input bam_in{buffer_stream, format_bam{}};

        auto it = bam_in.begin();
        for (unsigned i = 0; i < 3; ++i, ++it)
        {
            EXPECT_EQ(get<field::seq>(*it), seqs[i]);
            EXPECT_EQ(get<field::id>(*it), ids[i]);
            EXPECT_EQ(get<field::qual>(*it), quals[i]);
        }
    }
};

TEST_F(read_bam, tags)
{
    std::string input
    {
R"(ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*+,-./++-	FI:i:1
ID2	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE	AS:i:3
ID3 lala	0	*	0	0	*	*	0	0	ACGTTTA	!!!!!!!	TI:i:2
)"
    };
    do_read_test(input);
}

TEST_F(read_bam, no_qual)
{
   std::string input
   {
R"(ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	*
ID2	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE
ID3 lala	0	*	0	0	*	*	0	0	ACGTTTA	!!!!!!!
)"
    };

    quals[0] = { "!!!!!!!!!!!!!!!!!!"_phred42 }; // filled with 0xFF in bam if not present.
    do_read_test(input);
}

TEST_F(read_bam, qual_too_short)
{
    // This test removes a char from the quality but keeps the block size the same.
    std::stringstream istream{single_record_bam_raw.erase(89, 1)};

    sequence_file_input fin{istream, format_bam{}};
    EXPECT_THROW(fin.begin(), unexpected_end_of_input);
}

TEST_F(read_bam, qual_too_long)
{
    // This test adds a char to the quality but keeps the block size the same.
    single_record_bam_raw.push_back('\x0C');
    std::stringstream istream{single_record_bam_raw};

    sequence_file_input fin{istream, format_bam{}};

    // expect that some error is thrown here.
    auto parse_file = [] (auto & fin)
    {
        for (auto && record : fin)
            (void) record;
    };

    EXPECT_THROW(parse_file(fin), unexpected_end_of_input);
}

TEST_F(read_bam, no_seq)
{
    // This test removes the sequence but keeps the block size the same.
    single_record_bam_raw.erase(63, 9)[43] = '\x00';  // remove sequence from raw string.
    std::stringstream istream{single_record_bam_raw};

    sequence_file_input fin{istream, format_bam{}};

    // expect that some error is thrown here.
    auto parse_file = [] (auto & fin)
    {
        for (auto && record : fin)
            (void) record;
    };
    EXPECT_THROW(parse_file(fin), unexpected_end_of_input);
}

// ----------------------------------------------------------------------------
// writing
// ----------------------------------------------------------------------------
struct write : public ::testing::Test
{
    std::string single_record_bam_raw
    {
        // First comment is byte offset of the first byte.
        // Header - 23 byte
        /*00*/ '\x42', '\x41', '\x4D', '\x01',                                             // BAM1
        /*04*/ '\x0B', '\x00', '\x00', '\x00',                                             // header size (11)
        /*08*/ '\x40', '\x48', '\x44', '\x09', '\x56', '\x4E', '\x3A', '\x31', '\x2E', '\x36', '\x0A', // @HD VN:1.6
        /*19*/ '\x00', '\x00', '\x00', '\x00',                                             // n_ref (0)
        // Alignment record - 67 byte
        /*23*/ '\x3F', '\x00', '\x00', '\x00',                                             // block_size (63)
        /*27*/ '\xFF', '\xFF', '\xFF', '\xFF',                                             // refID (-1)
        /*31*/ '\xFF', '\xFF', '\xFF', '\xFF',                                             // pos (-1)
        /*35*/ '\x04',                                                                     // l_read_name (4)
        /*36*/ '\x00',                                                                     // mapq (0)
        /*37*/ '\x48', '\x12',                                                             // bin (4680)
        /*39*/ '\x00', '\x00',                                                             // n_cigar_bin (0)
        /*41*/ '\x00', '\x00',                                                             // flag (0)
        /*43*/ '\x12', '\x00', '\x00', '\x00',                                             // l_seq (18)
        /*47*/ '\xFF', '\xFF', '\xFF', '\xFF',                                             // next_refId (-1)
        /*51*/ '\xFF', '\xFF', '\xFF', '\xFF',                                             // next_pos (-1)
        /*55*/ '\x00', '\x00', '\x00', '\x00',                                             // tlen (0)
        /*59*/ '\x49', '\x44', '\x31', '\x00',                                             // read_name (ID1)
                                                                                           // cigar (*)
        /*63*/ '\x12', '\x48', '\x88', '\x88', '\x88', '\x88', '\x88', '\x88', '\x88',     // seq  (ACGTTTTTTTTTTTTTTT)
        /*72*/ '\x00', '\x02', '\x02', '\x03', '\x04', '\x05', '\x06', '\x07', '\x08',     // qual (!##$%&'()*+,-./++-)
        /*81*/ '\x09', '\x0A', '\x0B', '\x0C', '\x0D', '\x0E', '\x0A', '\x0A', '\x0C'
    };

    dna5_vector seq{"ACGTTTTTTTTTTTTTTT"_dna5};
    std::string id{"ID1"};
    std::vector<phred42> qual{"!##$%&'()*+,-./++-"_phred42};

    sequence_file_output_options options{};
    std::ostringstream ostream;

    void do_write_test()
    {
    	sequence_file_output fout{ostream, format_bam{}};
        EXPECT_NO_THROW(( fout.emplace_back(seq, id, qual) ));
        ostream.flush();
    }
};

TEST_F(write, no_id)
{
    single_record_bam_raw[23] = '\x3D'; // Update the record size.
    single_record_bam_raw[35] = '\x02'; // If id is empty bam writes "*\0".
    single_record_bam_raw[59] = '\x2A'; // '*'
    single_record_bam_raw[60] = '\x00'; // '\0'
    single_record_bam_raw.erase(61, 2); // Erase the remainign bytes from the standard name.

    id = ""; // zero the id.
    do_write_test();

    EXPECT_EQ(ostream.str(), single_record_bam_raw);
}

TEST_F(write, with_no_seq_and_no_qual)
{
    single_record_bam_raw[23] = '\x24'; // Update the record size.
    single_record_bam_raw[43] = '\x00'; // Set l_seq to 0.
    single_record_bam_raw.erase(63, 27); // Remove sequence and quality.

    seq = ""_dna5;
    qual = ""_phred42;
    do_write_test();

    EXPECT_EQ(ostream.str(), single_record_bam_raw);
}

TEST_F(write, with_no_seq_but_qual)
{
    single_record_bam_raw[23] = '\x24'; // Update the record size.
    single_record_bam_raw[43] = '\x00'; // Set l_seq to 0.
    single_record_bam_raw.erase(63, 9); // Remove sequence but not quality.

    seq = ""_dna5;
    sequence_file_output fout{ostream, format_bam{}};
    EXPECT_THROW(( fout.emplace_back(seq, id, qual) ), format_error);
}

// // No qualities given
TEST_F(write, with_seq_but_no_qual)
{
    single_record_bam_raw.replace(72, 18, 18, '\xFF'); // replace quality with \FF.

    qual = ""_phred42;
    do_write_test();

    EXPECT_EQ(ostream.str(), single_record_bam_raw);
}
