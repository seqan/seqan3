// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <sstream>
#include <string>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/test/pretty_printing.hpp>

#include "sequence_file_format_test_template.hpp"

using namespace seqan3;

template <>
struct sequence_file_read<format_genbank> : public sequence_file_data
{
    std::string standard_input
    {
R"(LOCUS ID1
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
LOCUS ID3 lala
DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
            cds.
ACCESSION   ID3
ORIGIN
        1 ACGTTTA
//
)"
    };

    std::string illegal_alphabet_character_input
    {
R"(LOCUS ID1
DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
            cds.
ACCESSION   ID1
ORIGIN
        1 ACGTTTT?TT TTTTTTTT
//
)"
    };

    std::string standard_output
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
LOCUS       ID3 lala                 7 bp
ORIGIN
        1 ACGTTTA
//
)"
    };

    std::string no_or_ill_formatted_id_input
    {
R"(LOCOS ID1    stuff
ORIGIN
        1 ACGTTTTTTT TTTTTTTT
//
)"
    };
};

// ---------------------------------------------------------------------------------------------------------------------
// parametrized tests
// ---------------------------------------------------------------------------------------------------------------------

INSTANTIATE_TYPED_TEST_SUITE_P(genbank, sequence_file_read, format_genbank, );
INSTANTIATE_TYPED_TEST_SUITE_P(genbank, sequence_file_write, format_genbank, );

// ----------------------------------------------------------------------------
// reading
// ----------------------------------------------------------------------------

struct read : public sequence_file_read<format_genbank>
{
    sequence_file_input_options<dna15, false> options{};

    void do_read_test(std::string const & input)
    {
        std::stringstream istream{input};

        sequence_file_input fin{istream, format_genbank{}};
        fin.options = options;

        auto it = fin.begin();
        for (unsigned i = 0; i < 3; ++i, ++it)
        {
            EXPECT_EQ(get<field::id>(*it), ids[i]);
            EXPECT_EQ(get<field::seq>(*it), seqs[i]);
        }
    }
};

TEST_F(read, complete_header)
{
    options.embl_genbank_complete_header = true;
    ids[0] = R"(LOCUS ID1
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
    ids[1] = R"(LOCUS ID2
DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
            cds.
ACCESSION   ID2
)";
    ids[2] = R"(LOCUS ID3 lala
DEFINITION  Homo sapiens mRNA for prepro cortistatin like peptide, complete
            cds.
ACCESSION   ID3
)";

    do_read_test(standard_input);
}

// ----------------------------------------------------------------------------
// writing
// ----------------------------------------------------------------------------

struct write : public sequence_file_data
{
    sequence_file_output_options options{};

    std::ostringstream ostream;

    void do_write_test()
    {
        sequence_file_output fout{ostream, format_genbank{}, fields<field::seq, field::id>{}};
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
