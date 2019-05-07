// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
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

using namespace seqan3;

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TEST(general, concepts)
{
    EXPECT_TRUE((SequenceFileInputFormat<sequence_file_format_embl>));
    EXPECT_TRUE((SequenceFileOutputFormat<sequence_file_format_embl>));
}

// ----------------------------------------------------------------------------
// reading
// ----------------------------------------------------------------------------

struct read : public ::testing::Test
{
    std::vector<std::string> expected_ids
    {
        { "ID1" },
        { "ID2" },
        { "ID3 lala" },
    };

    std::vector<dna5_vector> expected_seqs
    {
        { "ACGTTTTTTTTTTTTTTT"_dna5 },
        { "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"_dna5 },
        { "ACGTTTA"_dna5 },
    };

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

    sequence_file_format_embl format;

    sequence_file_input_options<dna5, false> options;

    std::string id;
    dna5_vector seq;

    void do_read_test(std::string const & input)
    {
        std::stringstream istream{input};

        for (unsigned i = 0; i < 3; ++i)
        {
            id.clear();
            seq.clear();

            EXPECT_NO_THROW(( format.read(istream, options, seq, id, std::ignore) ));
            EXPECT_EQ(id, expected_ids[i]);
            EXPECT_EQ(seq, expected_seqs[i]);
            EXPECT_TRUE((ranges::equal(seq, expected_seqs[i])));
            EXPECT_TRUE((ranges::equal(id, expected_ids[i])));
        }
    }
};

TEST_F(read, standard)
{
    do_read_test(input);
}

TEST_F(read, no_id)
{
    std::string input
    {
R"(IK ID1;  stuff
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

    std::stringstream istream{input};
    EXPECT_THROW(( format.read(istream, options, seq, id, std::ignore)), parse_error );
}

TEST_F(read, options_truncate_ids)
{
    options.truncate_ids = true;
    expected_ids[2] = "ID3"; // "lala" is stripped
    do_read_test(input);
}

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
    expected_ids[0] = "ID ID1;\tstuff\n";
    expected_ids[1] = "ID ID2;\n";
    expected_ids[2] = "ID ID3 lala;\nXX\nAC   AB000263;\nXX\n";
    do_read_test(input);
}

TEST_F(read, only_seq)
{
    std::stringstream istream{input};

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq.clear();

        format.read(istream, options, seq, std::ignore, std::ignore);

        EXPECT_TRUE((ranges::equal(seq, expected_seqs[i])));
    }
}

TEST_F(read, only_seq_multiple_lines_before)
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

    std::stringstream istream{input};

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq.clear();

        format.read(istream, options, seq, std::ignore, std::ignore);

        EXPECT_TRUE((ranges::equal(seq, expected_seqs[i])));
    }
}

TEST_F(read, only_id)
{
    std::stringstream istream{input};

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq.clear();

        format.read(istream, options, std::ignore, id, std::ignore);

        EXPECT_TRUE((ranges::equal(id, expected_ids[i])));
    }
}

TEST_F(read, seq_qual)
{
    std::stringstream istream{input};
    sequence_file_input_options<dna5, true> options2;

    std::vector<qualified<dna5, phred42>> seq_qual;

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq_qual.clear();

        format.read(istream, options2, seq_qual, id, seq_qual);

        EXPECT_TRUE((ranges::equal(id, expected_ids[i])));
        EXPECT_TRUE((ranges::equal(seq_qual | view::convert<dna5>, expected_seqs[i])));
    }
}

TEST_F(read, illegal_alphabet)
{
    std::string input
    {
        R"(ID ID1;	stuff
        SQ Sequence 1859 BP; 609 A; 314 C; 355 G; 581 T; 0 other;
          ARGTTTTTTT TTTTTTTT        18
        //)"
    };

    std::stringstream istream{input};
    EXPECT_THROW(( format.read(istream, options, seq, id, std::ignore)), parse_error );
}

TEST_F(read, from_stream_file)
{
    sequence_file_input fin{std::istringstream{input}, sequence_file_format_embl{}, fields<field::SEQ, field::ID>{}};

    size_t counter = 0;
    for (auto & [ seq, id ] : fin)
    {
        EXPECT_TRUE((std::ranges::equal(seq,  expected_seqs[counter])));
        EXPECT_TRUE((std::ranges::equal(id,  expected_ids[counter])));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

// ----------------------------------------------------------------------------
// writing
// ----------------------------------------------------------------------------

struct write : public ::testing::Test
{
    std::vector<dna5_vector> seqs
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

    sequence_file_format_embl format;

    sequence_file_output_options options;

    std::ostringstream ostream;

    void do_write_test()
    {
        for (unsigned i = 0; i < 3; ++i)
            EXPECT_NO_THROW(( format.write(ostream, options, seqs[i], ids[i], std::ignore) ));

        ostream.flush();
    }
};

TEST_F(write, arg_handling_id_missing)
{
    EXPECT_THROW( (format.write(ostream, options, seqs[0], std::ignore, std::ignore)),
                   std::logic_error );
}

TEST_F(write, arg_handling_id_empty)
{
    EXPECT_THROW( (format.write(ostream, options, seqs[0], std::string_view{""}, std::ignore)),
                   std::runtime_error );
}

TEST_F(write, arg_handling_seq_missing)
{
    EXPECT_THROW( (format.write(ostream, options, std::ignore, ids[0], std::ignore)),
                   std::logic_error );
}

TEST_F(write, arg_handling_seq_empty)
{
    EXPECT_THROW( (format.write(ostream, options, std::string_view{""}, ids[0], std::ignore)),
                   std::runtime_error );
}

TEST_F(write, default_options)
{

    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, seq_qual)
{
    auto convert_to_qualified = ranges::view::transform([] (auto const in)
    {
        return qualified<dna5, phred42>{} = in;
    });

    for (unsigned i = 0; i < 3; ++i)
        EXPECT_NO_THROW(( format.write(ostream,
                                       options,
                                       seqs[i] | convert_to_qualified,
                                       ids[i],
                                       seqs[i] | convert_to_qualified) ));

    ostream.flush();

    EXPECT_EQ(ostream.str(), comp);
}

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

TEST_F(write, from_stream_file)
{
    sequence_file_output fout{std::ostringstream{}, sequence_file_format_embl{}};

    for(int i = 0; i < 3; i++)
    {
        fout.emplace_back(seqs[i],ids[i]);
    }

    fout.get_stream().flush();

    EXPECT_EQ(reinterpret_cast<std::ostringstream&>(fout.get_stream()).str(), comp);
}
