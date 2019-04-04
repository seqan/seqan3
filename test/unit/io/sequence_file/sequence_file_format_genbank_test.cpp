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
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/io/sequence_file/format_genbank.hpp>
#include <seqan3/range/view/convert.hpp>

using namespace seqan3;

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TEST(general, concepts)
{
    EXPECT_TRUE((sequence_file_input_format_concept<sequence_file_format_genbank>));
    EXPECT_TRUE((sequence_file_output_format_concept<sequence_file_format_genbank>));
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
        { "ID3" },
    };

    std::vector<dna5_vector> expected_seqs
    {
        { "ACGTTTTTTTTTTTTTTT"_dna5 },
        { "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"_dna5 },
        { "ACGTTTA"_dna5 },
    };

    sequence_file_format_genbank format;

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
    std::string input
    {
        "LOCUS ID1\tstuff\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID1\n"
        "ORIGIN\n"
        "        1 ACGTTTTTTT TTTTTTTT\n"
        "//\n"
        "LOCUS ID2\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID2\n"
        "ORIGIN\n"
        "        1  ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT\n"
        "       61 TTTTTTTTTT TTTTTTTTTT TT\n"
        "//\n"
        "LOCUS ID3\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID3\n"
        "ORIGIN\n"
        "        1 ACGTTTA\n"
        "//\n"
    };

    do_read_test(input);
}

TEST_F(read, whitespace_in_seq)
{
    std::string input
    {
        "LOCUS ID1\tstuff\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID1\n"
        "ORIGIN\n"
        "        1 ACGTTTTTTT TTTTT\vTTT\n"
        "//\n"
        "LOCUS ID2\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID2\n"
        "ORIGIN\n"
        "        1  ACGTTTTTTT TTTTTT\nTTTT TTTTTTTTTT TTTTTTTTTT\t TTTTTTTTTT TTTTTTTTTT\n"
        "       61 TTTTTTTTTT TTTTTTTTTT TT\n"
        "//\n"
        "LOCUS ID3\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID3\n"
        "ORIGIN\n"
        "        1 ACGT\tTTA\n"
        "//\n"
    };
    do_read_test(input);
}

TEST_F(read, complete_header)
{
    std::string input
    {
        "LOCUS ID1\tstuff\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID1\n"
        "ORIGIN\n"
        "        1 ACGTTTTTTT TTTTTTTT\n"
        "//\n"
        "LOCUS ID2\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID2\n"
        "ORIGIN\n"
        "        1  ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT\n"
        "       61 TTTTTTTTTT TTTTTTTTTT TT\n"
        "//\n"
        "LOCUS ID3\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID3\n"
        "ORIGIN\n"
        "        1 ACGTTTA\n"
        "//\n"
    };

    options.complete_header = true;
    expected_ids[0] = std::string{"LOCUS ID1\tstuff\nDEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide,"
    " complete\n            cds.\nACCESSION   ID1\n"};
    expected_ids[1] = std::string{"LOCUS ID2\nDEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide,"
    " complete\n            cds.\nACCESSION   ID2\n"};
    expected_ids[2] = std::string{"LOCUS ID3\nDEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide,"
    " complete\n            cds.\nACCESSION   ID3\n"};
    do_read_test(input);
}

TEST_F(read, only_seq)
{
    std::string input
    {
        "LOCUS ID1\tstuff\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID1\n"
        "ORIGIN\n"
        "        1 ACGTTTTTTT TTTTTTTT\n"
        "//\n"
        "LOCUS ID2\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID2\n"
        "ORIGIN\n"
        "        1  ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT\n"
        "       61 TTTTTTTTTT TTTTTTTTTT TT\n"
        "//\n"
        "LOCUS ID3\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID3\n"
        "ORIGIN\n"
        "        1 ACGTTTA\n"
        "//\n"
    };

    std::stringstream istream{input};

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq.clear();

        EXPECT_NO_THROW( (format.read(istream, options, seq, std::ignore, std::ignore) ));

        EXPECT_TRUE((ranges::equal(seq, expected_seqs[i])));
    }
}

TEST_F(read, only_id)
{
    std::string input
    {
        "LOCUS ID1\tstuff\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID1\n"
        "ORIGIN\n"
        "        1 ACGTTTTTTT TTTTTTTT\n"
        "//\n"
        "LOCUS ID2\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID2\n"
        "ORIGIN\n"
        "        1  ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT\n"
        "       61 TTTTTTTTTT TTTTTTTTTT TT\n"
        "//\n"
        "LOCUS ID3\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID3\n"
        "ORIGIN\n"
        "        1 ACGTTTA\n"
        "//\n"
    };

    std::stringstream istream{input};

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq.clear();

        EXPECT_NO_THROW( (format.read(istream, options, std::ignore, id, std::ignore) ));

        EXPECT_TRUE((ranges::equal(id, expected_ids[i])));
    }
}

TEST_F(read, ignore_id)
{
    std::string input
    {
        "LOCUS ID1\tstuff\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID1\n"
        "ORIGIN\n"
        "        1 ACGTTTTTTT TTTTTTTT\n"
        "//\n"
        "LOCUS ID2\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID2\n"
        "ORIGIN\n"
        "        1  ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT\n"
        "       61 TTTTTTTTTT TTTTTTTTTT TT\n"
        "//\n"
        "LOCUS ID3\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID3\n"
        "ORIGIN\n"
        "        1 ACGTTTA\n"
        "//\n"
    };

    std::stringstream istream{input};

    for (unsigned i = 0; i < 3; ++i)
    {
        seq.clear();

        EXPECT_NO_THROW( (format.read(istream, options, seq, std::ignore, std::ignore) ));

        EXPECT_TRUE((ranges::equal(seq, expected_seqs[i])));
    }
}

TEST_F(read, no_locus)
{
    std::string input
    {
        "LOCOS ID1\tstuff\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID1\n"
        "ORIGIN\n"
        "        1 ACGTTTTTTT TTTTTTTT\n"
        "//\n"
    };

    std::stringstream istream{input};
    seq.clear();
    EXPECT_THROW( (format.read(istream, options, seq, std::ignore, std::ignore)), parse_error);

}

TEST_F(read, seq_qual)
{
    std::string input
    {
        "LOCUS ID1\tstuff\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID1\n"
        "ORIGIN\n"
        "        1 ACGTTTTTTT TTTTTTTT\n"
        "//\n"
        "LOCUS ID2\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID2\n"
        "ORIGIN\n"
        "        1  ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT\n"
        "       61 TTTTTTTTTT TTTTTTTTTT TT\n"
        "//\n"
        "LOCUS ID3\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID3\n"
        "ORIGIN\n"
        "        1 ACGTTTA\n"
        "//\n"
    };

    std::stringstream istream{input};
    sequence_file_input_options<dna5, true> options2;

    std::vector<qualified<dna5, phred42>> seq_qual;

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq_qual.clear();

        EXPECT_NO_THROW( (format.read(istream, options2, seq_qual, id, seq_qual) ));

        EXPECT_TRUE((ranges::equal(id, expected_ids[i])));
        EXPECT_TRUE((ranges::equal(seq_qual | view::convert<dna5>, expected_seqs[i])));
    }
}

TEST_F(read, illegal_alphabet)
{
    std::string input
    {
        "LOCUS ID1\tstuff\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID1\n"
        "ORIGIN\n"
        "        1 ACGTTTT?TT TTTTTTTT\n"
        "//\n"
    };

    std::stringstream istream{input};
    EXPECT_THROW(( format.read(istream, options, seq, id, std::ignore)), parse_error );
}

// ----------------------------------------------------------------------------
// writing
// ----------------------------------------------------------------------------

struct write : public ::testing::Test
{
    std::vector<dna5_vector> seqs
    {
        { "ACGTTTTTTTTTTTTTTT"_dna5 },
        { "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"_dna5 },
        { "ACGTTTA"_dna5 },
    };

    std::vector<std::string> ids
    {
        "ID1",
        "ID2",
        "ID3"
    };

    sequence_file_format_genbank format;

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
    std::string comp
    {
        "LOCUS       ID1                 18 bp\n"
        "ORIGIN\n"
        "        1 ACGTTTTTTT TTTTTTTT\n"
        "//\n"
        "LOCUS       ID2                 82 bp\n"
        "ORIGIN\n"
        "        1 ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT\n"
        "       61 TTTTTTTTTT TTTTTTTTTT TT\n"
        "//\n"
        "LOCUS       ID3                 7 bp\n"
        "ORIGIN\n"
        "        1 ACGTTTA\n"
        "//\n"
    };

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

    std::string comp
    {
        "LOCUS       ID1                 18 bp\n"
        "ORIGIN\n"
        "        1 ACGTTTTTTT TTTTTTTT\n"
        "//\n"
        "LOCUS       ID2                 82 bp\n"
        "ORIGIN\n"
        "        1 ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT\n"
        "       61 TTTTTTTTTT TTTTTTTTTT TT\n"
        "//\n"
        "LOCUS       ID3                 7 bp\n"
        "ORIGIN\n"
        "        1 ACGTTTA\n"
        "//\n"
    };

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, complete_header)
{
    std::string comp
    {
        "LOCUS       ID1                 18 bp\n"
        "DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete\n"
        "            cds.\n"
        "ACCESSION   ID1\n"
        "VERSION     ID1\n"
        "KEYWORDS    .\n"
        "SOURCE      .\n"
        "  ORGANISM  .\n"
        "ORIGIN\n"
        "        1 ACGTTTTTTT TTTTTTTT\n"
        "//\n"
        "LOCUS       ID2                 82 bp\n"
        "DEFINITION  ID2\n"
        "ACCESSION   ID2\n"
        "VERSION     ID2\n"
        "KEYWORDS    .\n"
        "SOURCE      .\n"
        "  ORGANISM  .\n"
        "ORIGIN\n"
        "        1 ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT\n"
        "       61 TTTTTTTTTT TTTTTTTTTT TT\n"
        "//\n"
        "LOCUS       ID3                 7 bp\n"
        "DEFINITION  ID3\n"
        "ACCESSION   ID3\n"
        "VERSION     ID3\n"
        "KEYWORDS    .\n"
        "SOURCE      .\n"
        "  ORGANISM  .\n"
        "ORIGIN\n"
        "        1 ACGTTTA\n"
        "//\n"
    };

    options.complete_header = true;
    ids[0] = std::string{"LOCUS       ID1                 18 bp\nDEFINITION  Homo sapiens mRNA for prepro cortistatin "
    "like peptide, complete\n            cds.\nACCESSION   ID1\nVERSION     ID1\nKEYWORDS    .\nSOURCE      .\n"
    "  ORGANISM  .\n"};
    ids[1] = std::string{"LOCUS       ID2                 82 bp\nDEFINITION  ID2\n"
    "ACCESSION   ID2\nVERSION     ID2\nKEYWORDS    .\nSOURCE      .\n  ORGANISM  .\n"};
    ids[2] = std::string{"LOCUS       ID3                 7 bp\nDEFINITION  ID3\n"
    "ACCESSION   ID3\nVERSION     ID3\nKEYWORDS    .\nSOURCE      .\n  ORGANISM  .\n"};
    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}
