// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/io/sequence_file/format_genbank.hpp>
#include <seqan3/test/pretty_printing.hpp>

using namespace seqan3;

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TEST(general, concepts)
{
    EXPECT_TRUE((SequenceFileInputFormat<sequence_file_format_genbank>));
    EXPECT_TRUE((SequenceFileOutputFormat<sequence_file_format_genbank>));
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

    std::string input
    {
R"(LOCUS ID1 stuff
DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
            cds.
ACCESSION   ID1
ACCESSION   U49845
VERSION     U49845.1  GI:1293613
KEYWORDS    .
SOURCE      Saccharomyces cerevisiae (baker's yeast)
  ORGANISM  Saccharomyces cerevisiae
            Eukaryota; Fungi; Ascomycota; Saccharomycotina; Saccharomycetes;
            Saccharomycetales; Saccharomycetaceae; Saccharomyces.
REFERENCE   1  (bases 1 to 5028)
FEATURES             Location/Qualifiers
     source          1..5028
ORIGIN
        1 ACGTTTTTTT TTTTTTTT
//
LOCUS ID2
DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
            cds.
ACCESSION   ID2
ORIGIN
        1  ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT
       61 TTTTTTTTTT TTTTTTTTTT TT
//
LOCUS ID3
DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
            cds.
ACCESSION   ID3
ORIGIN
        1 ACGTTTA
//
)"
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
    do_read_test(input);
}

TEST_F(read, complete_header)
{
    options.embl_genbank_complete_header = true;
    expected_ids[0] = R"(LOCUS ID1 stuff
DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
            cds.
ACCESSION   ID1
ACCESSION   U49845
VERSION     U49845.1  GI:1293613
KEYWORDS    .
SOURCE      Saccharomyces cerevisiae (baker's yeast)
  ORGANISM  Saccharomyces cerevisiae
            Eukaryota; Fungi; Ascomycota; Saccharomycotina; Saccharomycetes;
            Saccharomycetales; Saccharomycetaceae; Saccharomyces.
REFERENCE   1  (bases 1 to 5028)
FEATURES             Location/Qualifiers
     source          1..5028
)";
    expected_ids[1] = R"(LOCUS ID2
DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
            cds.
ACCESSION   ID2
)";
    expected_ids[2] = R"(LOCUS ID3
DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
            cds.
ACCESSION   ID3
)";
    do_read_test(input);
}

TEST_F(read, only_seq)
{
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
R"(LOCOS ID1	stuff
DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
            cds.
ACCESSION   ID1
ORIGIN
        1 ACGTTTTTTT TTTTTTTT
//
)"
    };

    std::stringstream istream{input};
    seq.clear();
    EXPECT_THROW( (format.read(istream, options, seq, std::ignore, std::ignore)), parse_error);

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

        EXPECT_NO_THROW( (format.read(istream, options2, seq_qual, id, seq_qual) ));

        EXPECT_TRUE((ranges::equal(id, expected_ids[i])));
        EXPECT_TRUE((ranges::equal(seq_qual | view::convert<dna5>, expected_seqs[i])));
    }
}

TEST_F(read, illegal_alphabet)
{
    std::string input
    {
R"(LOCUS ID1	stuff
DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
            cds.
ACCESSION   ID1
ORIGIN
        1 ACGTTTT?TT TTTTTTTT
//
)"
    };

    std::stringstream istream{input};
    EXPECT_THROW(( format.read(istream, options, seq, id, std::ignore)), parse_error );
}

TEST_F(read, from_stream_file)
{
    sequence_file_input fin{std::istringstream{input}, sequence_file_format_genbank{}, fields<field::SEQ, field::ID>{}};

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

    std::string comp
    {
R"(LOCUS       ID1                 18 bp
ORIGIN
        1 ACGTTTTTTT TTTTTTTT
//
LOCUS       ID2                 82 bp
ORIGIN
        1 ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT
       61 TTTTTTTTTT TTTTTTTTTT TT
//
LOCUS       ID3                 7 bp
ORIGIN
        1 ACGTTTA
//
)"
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
R"(LOCUS       ID1                 18 bp
DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
            cds.
ACCESSION   ID1
VERSION     ID1
KEYWORDS    .
SOURCE      .
  ORGANISM  .
ORIGIN
        1 ACGTTTTTTT TTTTTTTT
//
LOCUS       ID2                 82 bp
DEFINITION  ID2
ACCESSION   ID2
VERSION     ID2
KEYWORDS    .
SOURCE      .
  ORGANISM  .
ORIGIN
        1 ACGTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT TTTTTTTTTT
       61 TTTTTTTTTT TTTTTTTTTT TT
//
LOCUS       ID3                 7 bp
DEFINITION  ID3
ACCESSION   ID3
VERSION     ID3
KEYWORDS    .
SOURCE      .
  ORGANISM  .
ORIGIN
        1 ACGTTTA
//
)"
    };

    options.embl_genbank_complete_header = true;
    ids[0] = R"(LOCUS       ID1                 18 bp
DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
            cds.
ACCESSION   ID1
VERSION     ID1
KEYWORDS    .
SOURCE      .
  ORGANISM  .
)";
    ids[1] = R"(LOCUS       ID2                 82 bp
DEFINITION  ID2
ACCESSION   ID2
VERSION     ID2
KEYWORDS    .
SOURCE      .
  ORGANISM  .
)";
    ids[2] = R"(LOCUS       ID3                 7 bp
DEFINITION  ID3
ACCESSION   ID3
VERSION     ID3
KEYWORDS    .
SOURCE      .
  ORGANISM  .
)";
    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, from_stream_file)
{
    sequence_file_output fout{std::ostringstream{}, sequence_file_format_genbank{}};

    for(int i = 0; i < 3; i++)
    {
        fout.emplace_back(seqs[i],ids[i]);
    }

    fout.get_stream().flush();

    EXPECT_EQ(reinterpret_cast<std::ostringstream&>(fout.get_stream()).str(), comp);
}
