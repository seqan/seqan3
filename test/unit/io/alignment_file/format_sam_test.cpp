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
// #include <seqan3/io/alignment/alignment_file_in_format_concept.hpp>
#include <seqan3/io/alignment_file/OutputFormat.hpp>
#include <seqan3/io/alignment_file/format_sam.hpp>
#include <seqan3/range/view/convert.hpp>

using namespace seqan3;
using namespace seqan3::literal;

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
    std::vector<gapped<dna5>> ref_seq_gapped1 = {dna5::A, dna5::C, dna5::T, dna5::G};
    std::vector<gapped<dna5>> ref_seq_gapped2 = {dna5::A, dna5::C, dna5::T, dna5::G,
                                                dna5::A, dna5::T, dna5::C, dna5::G,
                                                dna5::A};
    std::vector<gapped<dna5>> ref_seq_gapped3 = {dna5::A, dna5::C, dna5::T, dna5::G,
                                                dna5::A, dna5::T, dna5::C, dna5::G};

    std::string ref_id = "ref";

    std::vector<unsigned> ref_offsets
    {
        0,
        1,
        2
    };

    std::vector<std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>>> alignments
    {
        {ref_seq_gapped1, std::vector<gapped<dna5>>{dna5::C, gap::GAP, dna5::G, dna5::T}},
        {ref_seq_gapped2, std::vector<gapped<dna5>>{dna5::G, dna5::G, dna5::G, dna5::C, dna5::T,
                                                    dna5::G, dna5::N, gap::GAP, dna5::A}},
        {ref_seq_gapped3, std::vector<gapped<dna5>>{dna5::G, gap::GAP, dna5::A, dna5::G,
                                                    dna5::T, dna5::A, gap::GAP, dna5::T}}
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
