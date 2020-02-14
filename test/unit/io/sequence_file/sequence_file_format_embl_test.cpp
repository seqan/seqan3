// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/io/sequence_file/format_embl.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/test/pretty_printing.hpp>

#include "sequence_file_format_test_template.hpp"

template <>
struct sequence_file_read<seqan3::format_embl> : public sequence_file_data
{
    std::string standard_input
    {
R"(ID ID1;	 stuff
SQ Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
  ACGTTTTTTT TTTTTTTT        18
//
ID ID2;
SQ Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
  ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT 60
TTTTTTTTTT TTTTTTTTTT TT        82
//
ID ID3 lala;
SQ Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
  ACGTTTA        7
//)"
    };

    std::string illegal_alphabet_character_input
    {
R"(ID ID1;	 stuff
SQ Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
  AXGTTTTTTT TTTTTTTT        18
//)"
    };

    std::string standard_output
    {
R"(ID ID1; 18 BP.
SQ Sequence 18 BP;
ACGTTTTTTT TTTTTTTT                                               18
//
ID ID2; 82 BP.
SQ Sequence 82 BP;
ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT 60
TTTTTTTTTT TTTTTTTTTT TT                                          82
//
ID ID3 lala; 7 BP.
SQ Sequence 7 BP;
ACGTTTA                                                           7
//
)"
    };

    std::string no_or_ill_formatted_id_input
    {
R"(IK ID1;   stuff
SQ Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
  ACGTTTTTTT TTTTTTTT        18
//)"
    };
};

// ---------------------------------------------------------------------------------------------------------------------
// parametrized tests
// ---------------------------------------------------------------------------------------------------------------------

INSTANTIATE_TYPED_TEST_SUITE_P(embl, sequence_file_read, seqan3::format_embl, );
INSTANTIATE_TYPED_TEST_SUITE_P(embl, sequence_file_write, seqan3::format_embl, );

// ----------------------------------------------------------------------------
// reading
// ----------------------------------------------------------------------------

struct read : public sequence_file_data
{
    std::string input
    {
R"(ID ID1;	stuff
SQ Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
  ACGTTTTTTT TTTTTTTT        18
//
ID ID2;
SQ Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
  ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT 60
TTTTTTTTTT TTTTTTTTTT TT        82
//
ID ID3 lala;
SQ Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
  ACGTTTA        7
//)"
    };

    seqan3::sequence_file_input_options<seqan3::dna15, false> options{};

    void do_read_test(std::string const & input)
    {
        std::stringstream istream{input};

        seqan3::sequence_file_input fin{istream, seqan3::format_embl{}, seqan3::fields<seqan3::field::id,
                                                                                       seqan3::field::seq>{}};
        fin.options = options;

        auto it = fin.begin();
        for (unsigned i = 0; i < 3; ++i, ++it)
        {
            EXPECT_EQ(seqan3::get<seqan3::field::id>(*it), ids[i]);
            EXPECT_EQ(seqan3::get<seqan3::field::seq>(*it), seqs[i]);
        }
    }
};

TEST_F(read, complete_header)
{
    std::string input
    {
R"(ID ID1;	stuff
SQ Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
  ACGTTTTTTT TTTTTTTT        18
//
ID ID2;
SQ Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
  ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT 60
TTTTTTTTTT TTTTTTTTTT TT        82
//
ID ID3 lala;
XX
AC   AB000263;
XX
SQ Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
  ACGTTTA        7
//)"
    };

    options.embl_genbank_complete_header = true;
    ids[0] = "ID ID1;\tstuff\n";
    ids[1] = "ID ID2;\n";
    ids[2] = "ID ID3 lala;\nXX\nAC   AB000263;\nXX\n";
    do_read_test(input);
}

TEST_F(read, multiple_lines_before_seq)
{
    std::string input
    {
R"(ID ID1;	stuff
XX
XX
XX
SQ Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
  ACGTTTTTTT TTTTTTTT        18
//
ID ID2;
SQ Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
  ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT 60
TTTTTTTTTT TTTTTTTTTT TT        82
//
ID ID3 lala;
SQ Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
  ACGTTTA        7
//)"
    };

    do_read_test(input);
}

// ----------------------------------------------------------------------------
// writing
// ----------------------------------------------------------------------------

struct write : public ::testing::Test
{
    std::vector<seqan3::dna5_vector> seqs
    {
        "ACGT"_dna5,
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN"_dna5,
        "GGAGTATAATATATATATATATAT"_dna5
    };

    std::vector<std::string> ids
    {
        "TEST 1",
        "Test2",
        "Test3"
    };

	std::string comp
    {
R"(ID TEST 1; 4 BP.
SQ Sequence 4 BP;
ACGT                                                              4
//
ID Test2; 91 BP.
SQ Sequence 91 BP;
AGGCTGNAGG CTGNAGGCTG NAGGCTGNAG GCTGNAGGCT GNAGGCTGNA GGCTGNAGGC 60
TGNAGGCTGN AGGCTGNAGG CTGNAGGCTG N                                91
//
ID Test3; 24 BP.
SQ Sequence 24 BP;
GGAGTATAAT ATATATATAT ATAT                                        24
//
)"
    };

    seqan3::sequence_file_output_options options{};

    std::ostringstream ostream;

    void do_write_test()
    {
        seqan3::sequence_file_output fout{ostream, seqan3::format_embl{}, seqan3::fields<seqan3::field::seq,
                                                                                         seqan3::field::id>{}};
        fout.options = options;

        for (unsigned i = 0; i < 3; ++i)
            EXPECT_NO_THROW((fout.emplace_back(seqs[i], ids[i])));

        ostream.flush();
    }
};

TEST_F(write, complete_header)
{
    std::string comp
    {
R"(ID TEST 1; 4 BP.
XX
SQ Sequence 4 BP;
ACGT                                                              4
//
ID Test2; 91 BP.
XX
SQ Sequence 91 BP;
AGGCTGNAGG CTGNAGGCTG NAGGCTGNAG GCTGNAGGCT GNAGGCTGNA GGCTGNAGGC 60
TGNAGGCTGN AGGCTGNAGG CTGNAGGCTG N                                91
//
ID Test3; 24 BP.
XX
SQ Sequence 24 BP;
GGAGTATAAT ATATATATAT ATAT                                        24
//
)"
    };
    options.embl_genbank_complete_header = true;
    ids[0] = std::string{"ID TEST 1; 4 BP.\nXX\n"};
    ids[1] = std::string{"ID Test2; 91 BP.\nXX\n"};
    ids[2] = std::string{"ID Test3; 24 BP.\nXX\n"};
    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}
