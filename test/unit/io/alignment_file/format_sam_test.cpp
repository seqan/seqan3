// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>
#include <range/v3/view/transform.hpp>

#include <seqan3/alphabet/quality/all.hpp>
// #include <seqan3/io/alignment/alignment_file_in_format_concept.hpp>
#include <seqan3/io/alignment_file/output_format_concept.hpp>
#include <seqan3/io/alignment_file/format_sam.hpp>
#include <seqan3/range/view/convert.hpp>

using namespace seqan3;

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TEST(general, concepts)
{
    // EXPECT_TRUE((alignment_file_in_format_concept<alignment_file_format_sam>));
    EXPECT_TRUE((AlignmentFileOutputFormat<alignment_file_format_sam>));
}

// ----------------------------------------------------------------------------
// writing
// ----------------------------------------------------------------------------

struct alignment_data : public ::testing::Test
{
    std::vector<dna5_vector> seqs
    {
        "ACGT"_dna5,
        "AGGCTGNAG"_dna5,
        "GGAGTATA"_dna5
    };

    std::vector<std::string> ids
    {
        "read1",
        "read2",
        "read3"
    };

    std::vector<std::vector<phred42>> quals
    {
        { "!##$"_phred42 },
        { "!##$&'()*"_phred42 },
        { "!!*+,-./"_phred42 },
    };

    std::vector<unsigned> offsets
    {
        1,
        0,
        1
    };

    dna5_vector ref_seq = "ACTGATCGAGAGGATCTAGAGGAGATCGTAGGAC"_dna5;
    std::vector<gapped<dna5>> ref_seq_gapped1 = {'A'_dna5, 'C'_dna5, 'T'_dna5, 'G'_dna5};
    std::vector<gapped<dna5>> ref_seq_gapped2 = {'A'_dna5, 'C'_dna5, 'T'_dna5, 'G'_dna5,
                                                'A'_dna5, 'T'_dna5, 'C'_dna5, 'G'_dna5,
                                                'A'_dna5};
    std::vector<gapped<dna5>> ref_seq_gapped3 = {'A'_dna5, 'C'_dna5, 'T'_dna5, 'G'_dna5,
                                                'A'_dna5, 'T'_dna5, 'C'_dna5, 'G'_dna5};

    std::string ref_id = "ref";

    std::vector<unsigned> ref_offsets
    {
        0,
        1,
        2
    };

    std::vector<std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>>> alignments
    {
        {ref_seq_gapped1, std::vector<gapped<dna5>>{'C'_dna5, gap{}, 'G'_dna5, 'T'_dna5}},
        {ref_seq_gapped2, std::vector<gapped<dna5>>{'G'_dna5, 'G'_dna5, 'G'_dna5, 'C'_dna5, 'T'_dna5,
                                                    'G'_dna5, 'N'_dna5, gap{}, 'A'_dna5}},
        {ref_seq_gapped3, std::vector<gapped<dna5>>{'G'_dna5, gap{}, 'A'_dna5, 'G'_dna5,
                                                    'T'_dna5, 'A'_dna5, gap{}, 'T'_dna5}}
    };

    std::vector<unsigned> flags
    {
        41,
        42,
        43
    };

    std::vector<unsigned> mapqs
    {
        61,
        62,
        63
    };

    std::vector<std::tuple<std::string, unsigned, unsigned>> mates
    {
        {"ref", 10, 300},
        {"ref", 10, 300},
        {"ref", 10, 300}
    };

    std::vector<sam_tag_dictionary> tag_dicts
    {
        sam_tag_dictionary{},
        sam_tag_dictionary{},
        sam_tag_dictionary{}
    };

};

TEST_F(alignment_data, default_options_all_members_specified)
{
    alignment_file_format_sam format;

    alignment_file_output_options options;

    std::ostringstream ostream;

    std::unique_ptr<alignment_file_header> header_ptr(nullptr);

    tag_dicts[0].get<"NM"_tag>() = 7;
    tag_dicts[0].get<"AS"_tag>() = 2;
    tag_dicts[1]["xy"_tag] = std::vector<uint16_t>{3,4,5};

    std::string comp
    {
        "read1\t41\tref\t1\t61\t1S1M1D2M\tref\t10\t300\tACGT\t!##$\tAS:i:2\tNM:i:7\n"
        "read2\t42\tref\t2\t62\t7M1D1M1S\tref\t10\t300\tAGGCTGNAG\t!##$&'()*\txy:B:S,3,4,5\n"
        "read3\t43\tref\t3\t63\t1S1M1D4M1D1M1S\tref\t10\t300\tGGAGTATA\t!!*+,-./\n"
    };

    for (unsigned i = 0; i < 3; ++i)
        /*EXPECT_NO_THROW(( */format.write(ostream, options, header_ptr,
                                       seqs[i], quals[i], ids[i], offsets[i], std::string{},
                                       ref_id, ref_offsets[i], alignments[i],
                                       flags[i], mapqs[i], mates[i],
                                       tag_dicts[i], 0u, 0u)/*))*/;

    ostream.flush();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(alignment_data, with_header)
{
    alignment_file_format_sam format;

    alignment_file_output_options options;

    std::ostringstream ostream;

    std::unique_ptr<alignment_file_header> header_ptr(new alignment_file_header);
    header_ptr->ref_dict["ref"] = {ref_seq.size(), ""};
    header_ptr->program_infos.push_back({"prog1", "cool_program", "", "", "", ""});
    header_ptr->comments.push_back("This is a comment.");

    tag_dicts[0].get<"NM"_tag>() = 7;
    tag_dicts[0].get<"AS"_tag>() = 2;
    tag_dicts[1]["xy"_tag] = std::vector<uint16_t>{3,4,5};

    std::string comp
    {
        "@HD\tVN:1.6\tSO:unknown\tGO:none\n"
        "@SQ\tSN:ref\tLN:34\n"
        "@PG\tID:prog1\tPN:cool_program\n"
        "@CO\tThis is a comment.\n"
        "read1\t41\tref\t1\t61\t1S1M1D2M\tref\t10\t300\tACGT\t!##$\tAS:i:2\tNM:i:7\n"
        "read2\t42\tref\t2\t62\t7M1D1M1S\tref\t10\t300\tAGGCTGNAG\t!##$&'()*\txy:B:S,3,4,5\n"
        "read3\t43\tref\t3\t63\t1S1M1D4M1D1M1S\tref\t10\t300\tGGAGTATA\t!!*+,-./\n"
    };

    for (unsigned i = 0; i < 3; ++i)
        EXPECT_NO_THROW(( format.write(ostream, options, header_ptr,
                                       seqs[i], quals[i], ids[i], offsets[i], std::string{},
                                       ref_id, ref_offsets[i], alignments[i],
                                       flags[i], mapqs[i], mates[i],
                                       tag_dicts[i], 0u, 0u)));

    ostream.flush();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(alignment_data, format_errors)
{
    alignment_file_format_sam format;
    alignment_file_output_options options;

    std::unique_ptr<alignment_file_header> header_ptr(new alignment_file_header);
    header_ptr->ref_dict["ref"] = {ref_seq.size(), ""};

    std::ostringstream ostream;

    // ensure that only a ref_id that is listed in the header is allowed
    EXPECT_THROW(( format.write(ostream, options, header_ptr,
                                   seqs[0], quals[0], ids[0], offsets[0], std::string{},
                                   "ref_id_that_does_not_exist", ref_offsets[0], alignments[0],
                                   flags[0], mapqs[0], mates[0],
                                   tag_dicts[0], 0u, 0u)),
                 format_error);
}
