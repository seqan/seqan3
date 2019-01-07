// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// ==========================================================================
#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>
#include <range/v3/view/transform.hpp>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/format_sam.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/range/view/convert.hpp>

using namespace seqan3;

inline std::vector<phred42> operator""_phred42(const char * s, std::size_t n)
{
    std::vector<phred42> r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}
// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------
TEST(general, concepts)
{
  EXPECT_TRUE((sequence_file_input_format_concept<sequence_file_format_sam>));
  EXPECT_TRUE((sequence_file_output_format_concept<sequence_file_format_sam>));
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
        "@ Comment\n"
        "@ Blablabla\n"
        "@ Bla\tbla\tbla\n"
        "ID1\t0\t*\t0\t0\t*\t*\t0\t0\tACGTTTTTTTTTTTTTTT\t!##$%&'()*+,-./++-\n"
        "ID2\t0\t*\t0\t0\t*\t*\t0\t0\tACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
        "TT\t!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
        "ID3 lala\t0\t*\t0\t0\t*\t*\t0\t0\tACGTTTA\t!!!!!!!\n"
    };
    do_read_test(input);
}

TEST_F(read, tags)
{
    std::string input
    {
      "@ Comment\n"
      "ID1\t0\tBABABA\t200\t0\t*\tBABABA\t0\t0\tACGTTTTTTTTTTTTTTT\t!##$%&'()*+,-./++-\tFI:i:1\n"
      "ID2\t0\t*\t0\t0\tBABA\t*\t30\t0\tACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
      "TTTTT\t!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\tAS:i:3\n"
      "ID3 lala\t0\t*\t0\t0\tHAHAHA+\t*\t0\t0\tACGTTTA\t!!!!!!!\tTI:i:2\n"
    };
    do_read_test(input);
}

TEST_F(read, whitespace_in_seq_qual)
{
    std::string input
    {
      "@ Comment\n"
      "ID1\t0\tBABABA\t200\t0\t*\tBABABA\t0\t0\tACGTTTTTTTTTTTTTTT\t!##$%&'()*+,-./++-\n"
      "ID2\t0\t*\t0\t0\tBABA\t*\t30\t0\tACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
      "TTTTTTT\t!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
      "ID3 lala\t0\t*\t0\t0\tHAHAHA+\t*\t0\t0\tACGTTTA\t!!!!!!!\n"
    };
    do_read_test(input);
}

TEST_F(read, mixed_issues)
{
  std::string input
  {
    "@ Comment\n"
    "ID1\t0\tBABABA\t200\t0\t*\tBABABA\t0\t0\tACGTTTTTTTTTTTTTTT\t!##$%&'()*+,-./++-\n"
    "ID2\t0\t*\t0\t0\tBABA\t*\t30\t0\tACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
    "TTTT\t!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
    "ID3 lala\t0\t*\t0\t0\tHAHAHA+\t*\t0\t0\tACGTTTA\t!!!!!!!\n"
    };
    do_read_test(input);
}

TEST_F(read, options_truncate_ids)
{
  std::string input
  {
      "@ Comment\n"
      "ID1 b\t0\t*\t0\t0\t*\t*\t0\t0\tACGTTTTTTTTTTTTTTT\t!##$%&'()*+,-./++-\n"
      "ID2 lala2\t0\t*\t0\t0\t*\t*\t0\t0\tACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
      "TTTTT\t!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
      "ID3 lala\t0\t*\t0\t0\t*\t*\t0\t0\tACGTTTA\t!!!!!!!\n"
    };

    options.truncate_ids = true;
    expected_ids[2] = "ID3"; // "lala" is stripped
    do_read_test(input);
}

TEST_F(read, no_qual)
{
   std::string input
   {
       "@ Comment\n"
       "ID1\t0\t*\t0\t0\t*\t*\t0\t0\tACGTTTTTTTTTTTTTTT\t*\n"
       "ID2\t0\t*\t0\t0\t*\t*\t0\t0\tACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
       "TTTTT\t!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
       "ID3 lala\t0\t*\t0\t0\t*\t*\t0\t0\tACGTTTA\t!!!!!!!\n"
    };

    expected_quals[0] = { ""_phred42 };
    do_read_test(input);
}

TEST_F(read, ignore_qual)
{
    std::string input
    {
        "@ Comment\n"
        "ID1\t0\t*\t0\t0\t*\t*\t0\t0\tACGTTTTTTTTTTTTTTT\t*\n"
        "ID2\t0\t*\t0\t0\t*\t*\t0\t0\tACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
        "TTTTT\t!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
        "ID3 lala\t0\t*\t0\t0\t*\t*\t0\t0\tACGTTTA\t!!!!!!!\n"
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
      "@ Comment\n"
      "ID 1\t0\tBABABA\t200\t0\t*\tBABABA\t0\t0\tACGTTTTTTTTTTTTTTT\t!##$%&'()-./++-\tFI:i:1\n"
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
      "@ Comment\n"
      "ID 1\t0\tBABABA\t200\t0\t*\tBABABA\t0\t0\tACGTTTTTTTTTTTTTTT\t!##$%&'()*+,-./++-+\tFI:i:1\n"
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
      "@ Comment\n"
      "ID 1\t0\tBABABA\t200\t0\t*\tBABABA\t0\t0\tACGTTTTTTTTTTTTTTT\t!##$%&'()*\n,-./++-\tFI:i:1\n"
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
      "@ Comment\n"
      "*\t0\tBABABA\t200\t0\t*\tBABABA\t0\t0\tACGTTTTTTTTTTTTTTT\t!##$%&'()*+,-./++-\tFI:i:1\n"
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
      "@ Comment\n"
      "ID 1\t0\tBABABA\t200\t0\t*\tBABABA\t0\t0\t*\t!##$%&'()*+,-./++-\tFI:i:1\n"
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
        "@ Comment\n"
        "ID1\t0\t*\t0\t0\t*\t*\t0\t0\tACGTTTTTTTTTTTTTTT\t!##$%&'()*+,-./++-\n"
        "ID2\t0\t*\t0\t0\t*\t*\t0\t0\tACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"
        "TT\t!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
        "ID3 lala\t0\t*\t0\t0\t*\t*\t0\t0\tACGTTTA\t!!!!!!!\n"
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
      "@ Comment\n"
      "ID 1\t0\tBABABA\t200\t0\t*\tBABABA\t0\t0\tACGTTTTT?TTTTTTTTT\t!##$%&'()*+,-./++-\tFI:i:1\n"
    };

    sequence_file_format_sam format;
    sequence_file_input_options<dna5, false> options;
    dna5_vector seq;
    std::stringstream istream{input};

    EXPECT_THROW((format.read(istream, options, seq, std::ignore, std::ignore)), parse_error );
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
      "TEST 1\t0\t*\t0\t0\t*\t*\t0\t0\tACGT\t!##$\n"
      "Test2\t0\t*\t0\t0\t*\t*\t0\t0\tAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGC"
      "TGNAGGCTGN\t!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE\n"
      "Test3\t0\t*\t0\t0\t*\t*\t0\t0\tGGAGTATAATATATATATATAT\t!!*+,-./+*+,-./+!!FF!!\n"
    };

    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, no_id)
{
    std::string comp
    {
      "*\t0\t*\t0\t0\t*\t*\t0\t0\tACGT\t!##$\n"
      "Test2\t0\t*\t0\t0\t*\t*\t0\t0\tAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGC"
      "TGNAGGCTGN\t!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE\n"
      "Test3\t0\t*\t0\t0\t*\t*\t0\t0\tGGAGTATAATATATATATATAT\t!!*+,-./+*+,-./+!!FF!!\n"
    };

    ids[0] = "";
    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, no_seq)
{
    std::string comp
    {
      "TEST 1\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"
      "Test2\t0\t*\t0\t0\t*\t*\t0\t0\tAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGC"
      "TGNAGGCTGN\t*\n"
      "Test3\t0\t*\t0\t0\t*\t*\t0\t0\tGGAGTATAATATATATATATAT\t*\n"
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
      "TEST 1\t0\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n"
      "Test2\t0\t*\t0\t0\t*\t*\t0\t0\tAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGC"
      "TGNAGGCTGN\t*\n"
      "Test3\t0\t*\t0\t0\t*\t*\t0\t0\tGGAGTATAATATATATATATAT\t*\n"
    };

    do_write_test_no_qual();

    EXPECT_EQ(ostream.str(), comp);
}
