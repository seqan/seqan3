// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/format_sam.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/test/pretty_printing.hpp>

using namespace seqan3;

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------
TEST(general, concepts)
{
  EXPECT_TRUE((SequenceFileInputFormat<sequence_file_format_sam>));
  EXPECT_TRUE((SequenceFileOutputFormat<sequence_file_format_sam>));
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
    std::vector<std::vector<phred42>> expected_quals
    {
        { "!##$%&'()*+,-./++-"_phred42 },
        { "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE"_phred42 },
        { "!!!!!!!"_phred42 },
    };

    sequence_file_format_sam format;
    sequence_file_input_options<dna5, false> options;
    std::string id;
    dna5_vector seq;
    std::vector<phred42> qual;
    void do_read_test(std::string const & input)
    {
        std::stringstream istream{input};
        for (unsigned i = 0; i < 3; ++i)
        {
            id.clear();
            seq.clear();
            qual.clear();
            EXPECT_NO_THROW(( format.read(istream, options, seq, id, qual) ));
            EXPECT_EQ(seq,expected_seqs[i]);
            EXPECT_EQ(id,expected_ids[i]);
            EXPECT_EQ(qual,expected_quals[i]);
        }
    }
};

TEST_F(read, standard)
{
    std::string input
    {
R"(@ Comment
@ Blablabla
@ Bla   bla bla
ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*+,-./++-
ID2	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE
ID3 lala	0	*	0	0	*	*	0	0	ACGTTTA	!!!!!!!
)"
    };
    do_read_test(input);
}

TEST_F(read, tags)
{
    std::string input
    {
R"(@ Comment
@ Blablabla
@ Bla   bla bla
ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*+,-./++-	FI:i:1
ID2	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE	AS:i:3
ID3 lala	0	*	0	0	*	*	0	0	ACGTTTA	!!!!!!!	TI:i:2
)"
    };
    do_read_test(input);
}

TEST_F(read, mixed_issues)
{
  std::string input
  {
R"(@ Comment
ID1	0	BABABA	200	0	*	BABABA	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*+,-./++-	FI:i:1
ID2	0	*	0	0	BABA	*	30	0	ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE
ID3 lala	0	*	0	0	HAHAHAHA+	*	0	0	ACGTTTA	!!!!!!!
)"
    };
    do_read_test(input);
}

TEST_F(read, options_truncate_ids)
{
  std::string input
  {
R"(@ Comment
ID1 b	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*+,-./++-
ID2 lala2	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE
ID3 lala	0	*	0	0	*	*	0	0	ACGTTTA	!!!!!!!
)"
    };

    options.truncate_ids = true;
    expected_ids[2] = "ID3"; // "lala" is stripped
    do_read_test(input);
}

TEST_F(read, seq_qual)
{

    std::string input
    {
R"(@ Comment
@ Blablabla
@ Bla   bla bla
ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*+,-./++-	FI:i:1
ID2	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE	AS:i:3
ID3 lala	0	*	0	0	*	*	0	0	ACGTTTA	!!!!!!!	TI:i:2
)"
    };

    std::stringstream istream{input};

    std::vector<qualified<dna5, phred42>> seq_qual;
    sequence_file_input_options<dna5, true> options2;

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq_qual.clear();

        format.read(istream, options2, seq_qual, id, seq_qual);

        EXPECT_TRUE((std::ranges::equal(id, expected_ids[i])));
        EXPECT_TRUE((std::ranges::equal(seq_qual | view::convert<dna5>, expected_seqs[i])));
        EXPECT_TRUE((std::ranges::equal(seq_qual | view::convert<phred42>, expected_quals[i])));
    }
}

TEST_F(read, no_qual)
{
   std::string input
   {
R"(@ Comment
ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	*
ID2	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE
ID3 lala	0	*	0	0	*	*	0	0	ACGTTTA	!!!!!!!
)"
    };

    expected_quals[0] = { ""_phred42 };
    do_read_test(input);
}

TEST_F(read, ignore_qual)
{
    std::string input
    {
R"(@ Comment
ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	*
ID2	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE
ID3 lala	0	*	0	0	*	*	0	0	ACGTTTA	!!!!!!!
)"
    };

    sequence_file_format_sam format;
    sequence_file_input_options<dna5, false> options;
    std::string id;
    dna5_vector seq;
    std::stringstream istream{input};

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq.clear();
        EXPECT_NO_THROW(( format.read(istream, options, seq, id, std::ignore) ));
        EXPECT_EQ(seq,expected_seqs[i]);
        EXPECT_EQ(id,expected_ids[i]);
    }
}

TEST_F(read, qual_too_short)
{
    std::string input
    {
R"(@ Comment
ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()-./++-
)"
    };

    sequence_file_format_sam format;
    sequence_file_input_options<dna5, false> options;
    std::string id;
    dna5_vector seq;
    std::vector<phred42> qual;
    std::stringstream istream{input};

    EXPECT_THROW((format.read(istream, options, seq, id, qual)), unexpected_end_of_input );
}

TEST_F(read, qual_too_long)
{
    std::string input
    {
R"(@ Comment
ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*+,-./++-+
)"
    };

    sequence_file_format_sam format;
    sequence_file_input_options<dna5, false> options;
    std::string id;
    dna5_vector seq;
    std::vector<phred42> qual;
    std::stringstream istream{input};

    EXPECT_THROW((format.read(istream, options, seq, id, qual)), unexpected_end_of_input );
}

TEST_F(read, wrong_qual3)
{
    std::string input
    {
R"(@ Comment
ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*\n,-./++-
)"
    };

    sequence_file_format_sam format;
    sequence_file_input_options<dna5, false> options;
    std::string id;
    dna5_vector seq;
    std::vector<phred42> qual;
    std::stringstream istream{input};

    EXPECT_THROW((format.read(istream, options, seq, id, qual)), unexpected_end_of_input );
}

TEST_F(read, no_id)
{
    std::string input
    {
R"(@ Comment
*	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*+,-./++-
)"
    };

    sequence_file_format_sam format;
    sequence_file_input_options<dna5, false> options;
    std::string id;
    std::stringstream istream{input};

    EXPECT_THROW((format.read(istream, options, std::ignore, id, std::ignore)), parse_error );
}

TEST_F(read, no_seq)
{
    std::string input
    {
R"(@ Comment
ID 1	0	*	0	0	*	*	0	0	*	!##$%&'()*+,-./++-
)"
    };

    sequence_file_format_sam format;
    sequence_file_input_options<dna5, false> options;
    dna5_vector seq;
    std::stringstream istream{input};

    EXPECT_THROW((format.read(istream, options, seq, std::ignore, std::ignore)), parse_error );
}

TEST_F(read, ignore_seq)
{
    std::string input
    {
R"(@ Comment
ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*+,-./++-
ID2	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE
ID3 lala	0	*	0	0	*	*	0	0	ACGTTTA	!!!!!!!
)"
    };

    sequence_file_format_sam format;
    sequence_file_input_options<dna5, false> options;
    std::string id;
    std::vector<phred42> qual;
    std::stringstream istream{input};

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        qual.clear();
        EXPECT_NO_THROW(( format.read(istream, options, std::ignore, id, qual) ));
        EXPECT_EQ(id,expected_ids[i]);
        EXPECT_EQ(qual,expected_quals[i]);
    }
}

TEST_F(read, wrong_seq)
{
    std::string input
    {
R"(@ Comment
ID 1	0	*	0	0	*	*	0	0	ACGTTTTT?TTTTTTTTT	!##$%&'()*+,-./++-
)"
    };

    sequence_file_format_sam format;
    sequence_file_input_options<dna5, false> options;
    dna5_vector seq;
    std::stringstream istream{input};

    EXPECT_THROW((format.read(istream, options, seq, std::ignore, std::ignore)), parse_error );
}

TEST_F(read, from_stream_file)
{
    std::string input
    {
R"(read1	41	ref	1	61	1S1M1D2M	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read 2	42	ref	2	62	7M1D1M1S	ref	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1D4M1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)"
    };

    std::vector<dna5_vector> seq_comp
    {
        "ACGT"_dna5,
        "AGGCTGNAG"_dna5,
        "GGAGTATA"_dna5
    };

    std::vector<std::string> id_comp
    {
        "read1",
        "read 2",
        "read3"
    };

    std::vector<std::vector<phred42>> qual_comp
    {
        { "!##$"_phred42 },
        { "!##$&'()*"_phred42 },
        { "!!*+,-./"_phred42 },
    };

    sequence_file_input fin{std::istringstream{input}, sequence_file_format_sam{}};

    size_t counter = 0;
    for (auto & [ seq, id, qual ] : fin)
    {
        EXPECT_TRUE((std::ranges::equal(seq,  seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(id,  id_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(qual,  qual_comp[counter])));

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
        "GGAGTATAATATATATATATAT"_dna5
    };
    std::vector<std::string> ids
    {
        "TEST 1",
        "Test2",
        "Test3"
    };
    std::vector<std::vector<phred42>> quals
    {
        { "!##$"_phred42 },
        { "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE"_phred42 },
        { "!!*+,-./+*+,-./+!!FF!!"_phred42 },
    };
    sequence_file_format_sam format;
    sequence_file_output_options options;
    std::ostringstream ostream;
    void do_write_test()
    {
        for (unsigned i = 0; i < 3; ++i)
            EXPECT_NO_THROW(( format.write(ostream, options, seqs[i], ids[i], quals[i]) ));
        ostream.flush();
    }
    void do_write_test_no_qual()
    {
        for (unsigned i = 0; i < 3; ++i)
            EXPECT_NO_THROW(( format.write(ostream, options, seqs[i], ids[i],""_phred42) ));
        ostream.flush();
    }
};
TEST_F(write, arg_handling_id_missing)
{
    EXPECT_NO_THROW( (format.write(ostream, options, seqs[0], std::ignore, quals[0])));
}
TEST_F(write, arg_handling_seq_missing)
{
    EXPECT_NO_THROW( (format.write(ostream, options, std::ignore, ids[0], quals[0])));
}
TEST_F(write, arg_handling_qual_missing)
{
    EXPECT_NO_THROW( (format.write(ostream, options, seqs[0], ids[0], std::ignore)));
}

TEST_F(write, default_options)
{
    std::string comp
    {
R"(TEST 1	0	*	0	0	*	*	0	0	ACGT	!##$
Test2	0	*	0	0	*	*	0	0	AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE
Test3	0	*	0	0	*	*	0	0	GGAGTATAATATATATATATAT	!!*+,-./+*+,-./+!!FF!!
)"
    };

    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, no_id)
{
    std::string comp
    {
R"(*	0	*	0	0	*	*	0	0	ACGT	!##$
Test2	0	*	0	0	*	*	0	0	AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE
Test3	0	*	0	0	*	*	0	0	GGAGTATAATATATATATATAT	!!*+,-./+*+,-./+!!FF!!
)"
    };

    ids[0] = "";
    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, no_seq)
{
    std::string comp
    {
R"(TEST 1	0	*	0	0	*	*	0	0	*	*
Test2	0	*	0	0	*	*	0	0	AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN	*
Test3	0	*	0	0	*	*	0	0	GGAGTATAATATATATATATAT	*
)"
    };

    seqs[0] = ""_dna5;
    do_write_test_no_qual();

    EXPECT_EQ(ostream.str(), comp);
}

// No qualities given
TEST_F(write, no_qual)
{
    std::string comp
    {
R"(TEST 1	0	*	0	0	*	*	0	0	ACGT	*
Test2	0	*	0	0	*	*	0	0	AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN	*
Test3	0	*	0	0	*	*	0	0	GGAGTATAATATATATATATAT	*
)"
    };

    do_write_test_no_qual();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, from_stream_file)
{
    sequence_file_output fout{std::ostringstream{}, sequence_file_format_sam{}};

    for(int i = 0; i < 3; i++)
    {
        fout.emplace_back(seqs[i],ids[i], quals[i]);
    }

    fout.get_stream().flush();

    std::string comp
    {
R"(TEST 1	0	*	0	0	*	*	0	0	ACGT	!##$
Test2	0	*	0	0	*	*	0	0	AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE
Test3	0	*	0	0	*	*	0	0	GGAGTATAATATATATATATAT	!!*+,-./+*+,-./+!!FF!!
)"
    };

    EXPECT_EQ(reinterpret_cast<std::ostringstream&>(fout.get_stream()).str(), comp);
}
