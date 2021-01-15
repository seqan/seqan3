// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>

#include "sam_file_format_test_template.hpp"

template <>
struct sam_file_read<seqan3::format_sam> : public sam_file_data
{
    // -----------------------------------------------------------------------------------------------------------------
    // formatted input
    // -----------------------------------------------------------------------------------------------------------------

	using stream_type = std::istringstream;

    std::string big_header_input{
R"(@HD	VN:1.6	SO:coordinate	SS:coordinate:queryname	GO:none
@PG	ID:qc	PN:quality_control	CL:qc -f file1	DS:trim reads with low qual	VN:1.0.0
@PG	ID:novoalign	PN:novoalign	VN:V3.02.07	CL:novoalign -d /path/hs37d5.ndx -f /path/file.fastq.gz	PP:qc
@SQ	SN:ref	LN:249250621
@SQ	SN:ref2	LN:243199373	AS:hs37d5
@RG	ID:U0a_A2_L1	PL:illumina	PU:1	LB:1	SM:NA12878
@RG	ID:U0a_A2_L2	PL:illumina	SM:NA12878	PU:1	LB:1
@CO	Tralalalalalala this is a comment
)"};

    std::string simple_three_reads_input{
R"(@HD	VN:1.6
@SQ	SN:ref	LN:34
read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read2	42	ref	2	62	1H7M1D1M1S2H	ref	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1P1M1I1M1I1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)"};

    std::string verbose_reads_input{
        "read1\t41\tref\t1\t61\t1S1M1D1M1I\t=\t10\t300\tACGT\t!##$\taa:A:c"
                                                                 "\tNM:i:-7"
                                                                 "\tAS:i:2"
                                                                 "\tff:f:3.1"
                                                                 "\tzz:Z:str"
                                                                 "\tCC:i:300"
                                                                 "\tcc:i:-300\n"
        "read2\t42\tref\t2\t62\t1H7M1D1M1S\tref\t10\t300\tAGGCTGNAG\t!##$&'()*\tbc:B:c,-3"
                                                                             "\tbC:B:C,3,200"
                                                                             "\tbs:B:s,-3,200,-300"
                                                                             "\tbS:B:S,300,40,500"
                                                                             "\tbi:B:i,-3,200,-66000"
                                                                             "\tbI:B:I,294967296"
                                                                             "\tbf:B:f,3.5,0.1,43.8"
                                                                             "\tbH:H:1AE301\n"
        "read3\t43\tref\t3\t63\t1S1M1P1M1I1M1I1D1M1S\tref\t10\t300\tGGAGTATA\t!!*+,-./\n"};

    std::string empty_input{"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:34\n*\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"};

    std::string empty_cigar{"read1\t41\tref\t1\t61\t*\tref\t10\t300\tACGT\t!##$\n"};

    std::string unknown_ref{"read1\t41\traf\t1\t61\t1S1M1D1M1I\t=\t10\t300\tACGT\t!##$\taa:A:c\tAS:i:2\tff:f:3.1\tzz:Z:str\n"};

    std::string unknown_ref_header{"@HD\tVN:1.6\n@SQ\tSN:ref\tLN:34\n*\t0\tunknown_ref\t1\t0\t4M\t*\t0\t0\tAAAA\t*\n"};

    std::string many_refs{[] ()
    {
        std::string result{"@HD\tVN:1.6\n"};
        for (size_t i = 0; i < 64; ++i)
            result += "@SQ\tSN:ref_" + std::to_string(i) + "\tLN:100\n";
        return result;
    }()};

    // -----------------------------------------------------------------------------------------------------------------
    // formatted output
    // -----------------------------------------------------------------------------------------------------------------

    std::string simple_three_reads_output{ // compared to simple_three_reads_input this has no hard clipping
R"(@HD	VN:1.6
@SQ	SN:ref	LN:34
read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read2	42	ref	2	62	7M1D1M1S	ref	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1P1M1I1M1I1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)"};

    std::string verbose_output{
R"(@HD	VN:1.6	SO:unknown	GO:none
@SQ	SN:ref	LN:34	AN:other_name
@RG	ID:group1	DS:more info
@PG	ID:prog1	PN:cool_program	CL:./prog1	PP:a	DS:b	VN:c
@CO	This is a comment.
read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	CC:i:300	NM:i:-7	aa:A:c	cc:i:-300	ff:f:3.1	zz:Z:str
read2	42	ref	2	62	7M1D1M1S	ref	10	300	AGGCTGNAG	!##$&'()*	bC:B:C,3,200	bI:B:I,294967296	bS:B:S,300,40,500	bc:B:c,-3	bf:B:f,3.5,0.1,43.8	bi:B:i,-3,200,-66000	bs:B:s,-3,200,-300
read3	43	ref	3	63	1S1M1P1M1I1M1I1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)"};

    std::string special_output{
R"(@HD	VN:1.6
@SQ	SN:ref	LN:34
read1	41	*	1	61	1S1M1D1M1I	*	0	0	ACGT	!##$
)"};

    std::string wrong_hexadecimal_tag{
R"(@SQ	SN:ref	LN:150
read1	41	ref	1	61	1S1M1D1M1I	=	10	300	ACGT	!##$	bH:H:1AE30
)"};
};

// ---------------------------------------------------------------------------------------------------------------------
// parametrized tests
// ---------------------------------------------------------------------------------------------------------------------

INSTANTIATE_TYPED_TEST_SUITE_P(sam, sam_file_read, seqan3::format_sam, );
INSTANTIATE_TYPED_TEST_SUITE_P(sam, sam_file_write, seqan3::format_sam, );

// ---------------------------------------------------------------------------------------------------------------------
// SAM specifics
// ---------------------------------------------------------------------------------------------------------------------

struct sam_format : public sam_file_data
{};

// since BAM uses the same read header function from SAM, it only needs to be tested once
TEST_F(sam_format, header_errors)
{
    {
        std::string header_str
        {
            "@HD\tVN:1.0\tTT:this is not a valid tag\n"
        };
        std::istringstream istream(header_str);
        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }
    {
        std::string header_str
        {
            "@HD\tVN:1.0\tSI:this is not a valid tag starting with S\n"
        };
        std::istringstream istream(header_str);
        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }
    {
        std::string header_str
        {
            "@HD\tVN:1.0\n"
            "@TT\tthis is not a valid tag\n"
        };
        std::istringstream istream(header_str);
        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }
    {
        std::string header_str
        {
            "@HD\tVN:1.0\n"
            "@PG\tID:prog\tTT:this is not a valid tag\n"
        };
        std::istringstream istream(header_str);
        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }
    {
        std::string header_str
        {
            "@HD\tVN:1.0\n"
            "@SQ\tSN:unknown_ref\tLN:0\n"
        };
        std::istringstream istream(header_str);
        seqan3::sam_file_input fin{istream, this->ref_ids, this->ref_sequences, seqan3::format_sam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }
    {
        std::string header_str
        {
            "@HD\tVN:1.0\n"
            "@SQ\tSN:ref\tLN:0\n" /*wrong length*/
        };
        std::istringstream istream(header_str);
        seqan3::sam_file_input fin{istream, this->ref_ids, this->ref_sequences, seqan3::format_sam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }
}

TEST_F(sam_format, no_hd_line_in_header)
{
    // the header line (@HD) is optional
    std::istringstream istream{std::string{"@SQ\tSN:ref\tLN:34\nread1\t41\tref\t1\t61\t*\tref\t10\t300\tACGT\t!##$\n"}};
    seqan3::sam_file_input fin{istream, seqan3::format_sam{}, seqan3::fields<seqan3::field::id>{}};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    EXPECT_EQ(seqan3::get<seqan3::field::id>(*fin.begin()), std::string{"read1"});
#pragma GCC diagnostic pop

    EXPECT_EQ((*fin.begin()).id(), std::string{"read1"});
}

TEST_F(sam_format, windows_file)
{
    std::istringstream istream(std::string("read1\t41\tref\t1\t61\t*\tref\t10\t300\tACGT\t!##$\r\n"));
    seqan3::sam_file_input fin{istream, seqan3::format_sam{}, seqan3::fields<seqan3::field::id>{}};

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    EXPECT_EQ(seqan3::get<seqan3::field::id>(*fin.begin()), std::string{"read1"});
#pragma GCC diagnostic pop

    EXPECT_EQ((*fin.begin()).id(), std::string{"read1"});
}

TEST_F(sam_format, format_error_illegal_character_in_seq)
{
    std::istringstream istream(std::string("*\t0\t*\t0\t0\t*\t*\t0\t0\tAC!T\t*\n"));

    seqan3::sam_file_input fin{istream, seqan3::format_sam{}};
    EXPECT_THROW(fin.begin(), seqan3::parse_error);
}

TEST_F(sam_format, format_error_invalid_arithmetic_value)
{
    // invalid value
    std::istringstream istream(std::string("*\t0\t*\t1abc\t0\t*\t*\t0\t0\t*\t*\n"));
    {
        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }
    // overflow error
    {
        istream = std::istringstream(std::string("*\t0\t*\t2147483650\t0\t*\t*\t0\t0\t*\t*\n"));
        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }
    // negative value as ref_offset
    {
        istream = std::istringstream(std::string("*\t0\t*\t-3\t0\t*\t*\t0\t0\t*\t*\n"));
        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }

    // negative value as mate mapping position
    {
        istream = std::istringstream(std::string("*\t0\t*\t0\t0\t*\t*\t-3\t0\t*\t*\n"));
        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }
}

TEST_F(sam_format, format_error_invalid_cigar)
{
    // unkown operation
    std::istringstream istream(std::string("*\t0\t*\t0\t0\t5Z\t*\t0\t0\t*\t*\n"));
    {
        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }
    // negative number as operation count
    {
        istream = std::istringstream(std::string("*\t0\t*\t0\t0\t-5M\t*\t0\t0\t*\t*\n"));
        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }

    {
        istream = std::istringstream(std::string("*\t0\t*\t0\t0\t3S4M1I-5M2D2M\t*\t0\t0\t*\t*\n"));
        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }
}

TEST_F(sam_format, format_error_invalid_sam_tag_format)
{
    // type identifier is wrong
    std::istringstream istream(std::string("*\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*\tNM:X:3\n"));
    {
        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }
    // Array subtype identifier is wrong
    {
        istream = std::istringstream(std::string("*\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*\tNM:B:x3,4\n"));
        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};
        EXPECT_THROW(fin.begin(), seqan3::format_error);
    }
}

TEST_F(sam_format, short_cigar_string_with_softclipping)
{
    using seqan3::operator""_dna5;

	// The member function transfer_soft_clipping_to needs to work on 2 element cigar strings
    {
        std::istringstream istream("id	16	ref	0	255	10M5S	*	0	0	AGAGGGGGATAACCA	*\n");
        seqan3::sam_file_input fin{istream,
                                   ref_ids,
                                   ref_sequences,
                                   seqan3::format_sam{},
                                   seqan3::fields<seqan3::field::alignment>{}};
        EXPECT_RANGE_EQ(std::get<1>(std::get<0>(*fin.begin())), "AGAGGGGGAT"_dna5);
    }

    {
        std::istringstream istream("id	16	ref	0	255	5S10M	*	0	0	AGAGGGGGATAACCA	*\n");
        seqan3::sam_file_input fin{istream,
                                   ref_ids,
                                   ref_sequences,
                                   seqan3::format_sam{},
                                   seqan3::fields<seqan3::field::alignment>{}};
        EXPECT_RANGE_EQ(std::get<1>(std::get<0>(*fin.begin())), "GGGATAACCA"_dna5);
    }
}

TEST_F(sam_format, write_different_header)
{
    std::ostringstream ostream;

    auto write_header = [&] ()
    {
        seqan3::sam_file_output fout{ostream, seqan3::format_sam{}, seqan3::fields<seqan3::field::header_ptr,
                                                                                   seqan3::field::ref_id,
                                                                                   seqan3::field::ref_offset>{}};
        ASSERT_NO_THROW(fout.emplace_back(&header, this->ref_id, 0));
    };

    header.sorting = "unsorted";
    header.grouping = "query";

    write_header();
    ostream.flush();
    EXPECT_EQ(ostream.str(),
              "@HD\tVN:1.6\tSO:unsorted\tGO:query\n@SQ\tSN:ref\tLN:34\n*\t0\tref\t1\t0\t*\t*\t0\t0\t*\t*\n");

    ostream = std::ostringstream{};
    header.sorting = "queryname";
    header.grouping = "reference";

    write_header();
    ostream.flush();
    EXPECT_EQ(ostream.str(),
              "@HD\tVN:1.6\tSO:queryname\tGO:reference\n@SQ\tSN:ref\tLN:34\n*\t0\tref\t1\t0\t*\t*\t0\t0\t*\t*\n");

    ostream = std::ostringstream{};
    header.sorting = "coordinate";
    header.subsorting = "query";

    write_header();
    ostream.flush();
    EXPECT_EQ(ostream.str(),
              "@HD\tVN:1.6\tSO:coordinate\tSS:query\tGO:reference\n@SQ\tSN:ref\tLN:34\n*\t0\tref\t1\t0\t*\t*\t0\t0\t*\t*\n");
}

TEST_F(sam_format, issue2195)
{ // see issue https://github.com/seqan/seqan3/issues/2195
    {
        std::istringstream istream
        {
            "*r1\t4\t1\t10\t0\t5M\t=\t136097\t-121\tACTGA\t*9<9;\tNM:i:1\tMQ:i:0\n"
        };
        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};

        using seqan3::operator""_phred42;
        std::vector<seqan3::phred42> expected_quality = "*9<9;"_phred42;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
        EXPECT_RANGE_EQ(seqan3::get<seqan3::field::id>(*fin.begin()), std::string{"*r1"});
        EXPECT_RANGE_EQ(seqan3::get<seqan3::field::qual>(*fin.begin()), expected_quality);
#pragma GCC diagnostic pop

        EXPECT_RANGE_EQ((*fin.begin()).id(), std::string{"*r1"});
        EXPECT_RANGE_EQ((*fin.begin()).base_qualities(), expected_quality);
    }

    {
        std::istringstream istream
        {
            "*\t4\t1\t10\t0\t2M\t=\t136097\t-121\tAC\t*1\tNM:i:1\tMQ:i:0\n"
        };
        seqan3::sam_file_input fin{istream, seqan3::format_sam{}};

        using seqan3::operator""_phred42;
        std::vector<seqan3::phred42> expected_quality = "*1"_phred42;

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
        EXPECT_RANGE_EQ(seqan3::get<seqan3::field::id>(*fin.begin()), std::string{""});
        EXPECT_RANGE_EQ(seqan3::get<seqan3::field::qual>(*fin.begin()), expected_quality);
#pragma GCC diagnostic pop

        EXPECT_RANGE_EQ((*fin.begin()).id(), std::string{""});
        EXPECT_RANGE_EQ((*fin.begin()).base_qualities(), expected_quality);
    }
}
