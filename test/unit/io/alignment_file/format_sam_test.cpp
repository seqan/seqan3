// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/alignment_file/input_format_concept.hpp>
#include <seqan3/io/alignment_file/output_format_concept.hpp>
#include <seqan3/io/alignment_file/format_sam.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/test/pretty_printing.hpp>

using namespace seqan3;

// global variables for reuse
alignment_file_input_options<dna5> input_options;
alignment_file_output_options output_options;

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TEST(general, concepts)
{
    EXPECT_TRUE((AlignmentFileOutputFormat<alignment_file_format_sam>));
    EXPECT_TRUE((AlignmentFileInputFormat<alignment_file_format_sam>));
}

struct reading_sam : public ::testing::Test
{
    reading_sam()
    {
        ref_sequences = std::vector<dna5_vector>{ref_seq};
        ref_ids = std::vector<std::string>{ref_id};
        header = alignment_file_header{ref_ids};
        header.ref_id_info.emplace_back(ref_seq.size(), "");
        header.ref_dict[header.ref_ids()[0]] = 0; // set up header which is otherwise done on file level
    }

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

    std::vector<int32_t> offsets
    {
        1,
        0,
        1
    };

    dna5_vector ref_seq = "ACTGATCGAGAGGATCTAGAGGAGATCGTAGGAC"_dna5;

    std::vector<gapped<dna5>> ref_seq_gapped1 = {'A'_dna5, 'C'_dna5, 'T'_dna5, gap{}};
    std::vector<gapped<dna5>> ref_seq_gapped2 = {'C'_dna5, 'T'_dna5, 'G'_dna5, 'A'_dna5,
                                                 'T'_dna5, 'C'_dna5, 'G'_dna5, 'A'_dna5, 'G'_dna5};
    std::vector<gapped<dna5>> ref_seq_gapped3 = {'T'_dna5, 'G'_dna5, 'A'_dna5, gap{},
                                                 'T'_dna5, gap{}, 'C'_dna5, 'G'_dna5,};

    std::string ref_id = "ref";

    std::vector<int32_t> ref_offsets
    {
        0,
        1,
        2
    };

    std::vector<std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>>> alignments
    {
        {ref_seq_gapped1, std::vector<gapped<dna5>>{'C'_dna5, gap{}, 'G'_dna5, 'T'_dna5}},
        {ref_seq_gapped2, std::vector<gapped<dna5>>{'A'_dna5, 'G'_dna5, 'G'_dna5, 'C'_dna5, 'T'_dna5,
                                                    'G'_dna5, 'N'_dna5, gap{}, 'A'_dna5}},
        {ref_seq_gapped3, std::vector<gapped<dna5>>{'G'_dna5, gap{}, 'A'_dna5, 'G'_dna5,
                                                    'T'_dna5, 'A'_dna5, gap{}, 'T'_dna5}}
    };

    std::vector<uint16_t> flags
    {
        41u,
        42u,
        43u
    };

    std::vector<uint8_t> mapqs
    {
        61u,
        62u,
        63u
    };

    std::vector<std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>> mates
    {
        {0, 9, 300},
        {0, 9, 300},
        {0, 9, 300}
    };

    std::vector<sam_tag_dictionary> tag_dicts
    {
        sam_tag_dictionary{},
        sam_tag_dictionary{},
        sam_tag_dictionary{}
    };

    std::vector<dna5_vector> ref_sequences{};
    std::vector<std::string> ref_ids{};
    alignment_file_header<std::vector<std::string>> header{};
};

// ----------------------------------------------------------------------------
// reading
// ----------------------------------------------------------------------------

struct read_header : public reading_sam
{};

TEST_F(read_header, sucess)
{
    alignment_file_format_sam format;

    std::string header_str =
R"(@HD	VN:1.0	SO:coordinate	SS:coordinate:queryname	GO:none
@PG	ID:qc	PN:quality_control	CL:qc -f file1	DS:trim reads with low qual	VN:1.0.0
@PG	ID:novoalign	PN:novoalign	VN:V3.02.07	CL:novoalign -d /path/hs37d5.ndx -f /path/file.fastq.gz	PP:qc
@SQ	SN:ref	LN:249250621
@SQ	SN:ref2	LN:243199373	AS:hs37d5
@RG	ID:U0a_A2_L1	PL:illumina	PU:1	LB:1	SM:NA12878
@RG	ID:U0a_A2_L2	PL:illumina	SM:NA12878	PU:1	LB:1
@CO	Tralalalalalala this is a comment
)";

    std::istringstream istream(header_str);

    alignment_file_header header{};

    ASSERT_NO_THROW(format.read(istream, input_options, std::ignore, header,  std::ignore, std::ignore, std::ignore,
                                std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                std::ignore, std::ignore, std::ignore, std::ignore, std::ignore));

    EXPECT_EQ(header.format_version, "1.0");
    EXPECT_EQ(header.sorting, "coordinate");
    EXPECT_EQ(header.subsorting, "coordinate:queryname");
    EXPECT_EQ(header.grouping, "none");

    EXPECT_EQ(header.program_infos[0].id, "qc");
    EXPECT_EQ(header.program_infos[0].name, "quality_control");
    EXPECT_EQ(header.program_infos[0].version, "1.0.0");
    EXPECT_EQ(header.program_infos[0].description, "trim reads with low qual");
    EXPECT_EQ(header.program_infos[0].previous, "");
    EXPECT_EQ(header.program_infos[0].command_line_call, "qc -f file1");
    EXPECT_EQ(header.program_infos[1].id, "novoalign");
    EXPECT_EQ(header.program_infos[1].name, "novoalign");
    EXPECT_EQ(header.program_infos[1].version, "V3.02.07");
    EXPECT_EQ(header.program_infos[1].description, "");
    EXPECT_EQ(header.program_infos[1].previous, "qc");
    EXPECT_EQ(header.program_infos[1].command_line_call, "novoalign -d /path/hs37d5.ndx -f /path/file.fastq.gz");

    std::string id1{"ref"};
    std::string id2{"ref2"};

    EXPECT_EQ(header.ref_id_info[(header.ref_dict[id1])], (std::tuple<uint32_t, std::string>{249250621u, ""}));
    EXPECT_EQ(header.ref_id_info[(header.ref_dict[id2])], (std::tuple<uint32_t, std::string>{243199373u, "AS:hs37d5"}));

    EXPECT_EQ(header.read_groups[0],
              (std::pair<std::string, std::string>{"U0a_A2_L1", "PL:illumina\tPU:1\tLB:1\tSM:NA12878"}));
    EXPECT_EQ(header.read_groups[1],
              (std::pair<std::string, std::string>{"U0a_A2_L2", "PL:illumina\tSM:NA12878\tPU:1\tLB:1"}));

    EXPECT_EQ(header.comments[0], "Tralalalalalala this is a comment");
}

TEST_F(read_header, errors)
{
    alignment_file_format_sam format;

    {
        std::string header_str
        {
            "@HD\tVN:1.0\tTT:this is not a valid tag\n"
        };
        std::istringstream istream(header_str);
        EXPECT_THROW(format.read(istream, input_options, std::ignore, header,  std::ignore, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                     format_error);
    }
    {
        std::string header_str
        {
            "@HD\tVN:1.0\tSI:this is not a valid tag starting with S\n"
        };
        std::istringstream istream(header_str);
        EXPECT_THROW(format.read(istream, input_options, std::ignore, header,  std::ignore, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                     format_error);
    }
    {
        std::string header_str
        {
            "@HD\tVN:1.0\n"
            "@TT\tthis is not a valid tag\n"
        };
        std::istringstream istream(header_str);
        EXPECT_THROW(format.read(istream, input_options, std::ignore, header,  std::ignore, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                     format_error);
    }
    {
        std::string header_str
        {
            "@HD\tVN:1.0\n"
            "@PG\tID:prog\tTT:this is not a valid tag\n"
        };
        std::istringstream istream(header_str);
        EXPECT_THROW(format.read(istream, input_options, std::ignore, header, std::ignore, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                     format_error);
    }
    {
        std::string header_str
        {
            "@HD\tVN:1.0\n"
            "@SQ\tSN:unknown_ref\tLN:0\n"
        };
        std::istringstream istream(header_str);
        EXPECT_THROW(format.read(istream, input_options, ref_sequences, header, std::ignore, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                     format_error);
    }
    {
        std::string header_str
        {
            "@HD\tVN:1.0\n"
            "@SQ\tSN:ref\tLN:0\n" /*wrong length*/
        };
        std::istringstream istream(header_str);
        EXPECT_THROW(format.read(istream, input_options, ref_sequences, header, std::ignore, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                     format_error);
    }
}

TEST_F(reading_sam, read_in_all_data)
{
    alignment_file_format_sam format;

    std::string file_in_str
    {
        "read1\t41\tref\t1\t61\t1S1M1D1M1I\t=\t10\t300\tACGT\t!##$\taa:A:c\tAS:i:2\tff:f:3.1\tzz:Z:str\n"
        "read2\t42\tref\t2\t62\t1H7M1D1M1S\tref\t10\t300\tAGGCTGNAG\t!##$&'()*\tbc:B:c,-3"
                                                                             "\tbC:B:C,3,200"
                                                                             "\tbs:B:s,-3,200,-300"
                                                                             "\tbS:B:S,300,40,500"
                                                                             "\tbi:B:i,-3,200,-66000"
                                                                             "\tbI:B:I,294967296"
                                                                             "\tbf:B:f,3.5,0.1,43.8\n"
        "read3\t43\tref\t3\t63\t1S1M1D1M1I1M1I1D1M1S\tref\t10\t300\tGGAGTATA\t!!*+,-./\n"
    };

    tag_dicts[0].get<"AS"_tag>() = 2;
    tag_dicts[0]["aa"_tag] = 'c';
    tag_dicts[0]["ff"_tag] = 3.1f;
    tag_dicts[0]["zz"_tag] = "str";
    tag_dicts[1]["bc"_tag] = std::vector<int8_t>{-3};
    tag_dicts[1]["bC"_tag] = std::vector<uint8_t>{3u, 200u};
    tag_dicts[1]["bs"_tag] = std::vector<int16_t>{-3, 200, -300};
    tag_dicts[1]["bS"_tag] = std::vector<uint16_t>{300u, 40u, 500u};
    tag_dicts[1]["bi"_tag] = std::vector<int32_t>{-3, 200, -66000};
    tag_dicts[1]["bI"_tag] = std::vector<uint32_t>{294967296u};
    tag_dicts[1]["bf"_tag] = std::vector<float>{3.5f, 0.1f, 43.8f};

    dna5_vector seq;
    std::string id;
    std::vector<phred42> qual;
    int32_t offset;
    std::optional<int32_t> ref_id_in;
    std::optional<int32_t> ref_offset;
    std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>> alignment;
    uint16_t flag;
    uint8_t mapq;
    std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> mate;
    sam_tag_dictionary tag_dict;

    std::istringstream istream(file_in_str);

    for (size_t i = 0; i < 3; ++i)
    {
        ASSERT_NO_THROW(format.read(istream, input_options, ref_sequences, header, seq, qual, id, offset, std::ignore,
                                    ref_id_in, ref_offset, alignment, flag, mapq, mate, tag_dict, std::ignore,
                                    std::ignore));

        EXPECT_EQ(seq, seqs[i]);
        EXPECT_EQ(id, ids[i]);
        EXPECT_EQ(qual, quals[i]);
        EXPECT_EQ(offset, offsets[i]);
        EXPECT_EQ(ref_id_in, 0);
        EXPECT_EQ(*ref_offset, ref_offsets[i]);
        EXPECT_EQ(get<0>(alignment), get<0>(alignments[i]));
        EXPECT_EQ(get<1>(alignment), get<1>(alignments[i]));
        EXPECT_EQ(flag, flags[i]);
        EXPECT_EQ(mapq, mapqs[i]);
        EXPECT_EQ(mate, mates[i]);
        EXPECT_EQ(tag_dict, tag_dicts[i]);

        seq.clear();
        id.clear();
        qual.clear();
        offset = 0;
        ref_id_in = 0;
        ref_offset = 0;
        alignment = std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>>{};
        flag = 0;
        mapq = 0;
        mate = std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>{};
        tag_dict.clear();
    }
}

TEST_F(reading_sam, read_in_all_but_empty_data)
{
    alignment_file_format_sam format;

    std::istringstream istream(std::string("*\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"));

    dna5_vector seq;
    std::string id;
    std::vector<phred42> qual;
    int32_t offset;
    std::optional<int32_t> ref_id_in;
    std::optional<int32_t> ref_offset;
    std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>> alignment;
    uint16_t flag;
    uint8_t mapq;
    std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> mate;
    sam_tag_dictionary tag_dict;

    ASSERT_NO_THROW(format.read(istream, input_options, ref_sequences, header, seq, qual, id, offset, std::ignore,
                                ref_id_in, ref_offset, alignment, flag, mapq, mate, tag_dict, std::ignore, std::ignore));

    EXPECT_TRUE(seq.empty());
    EXPECT_TRUE(id.empty());
    EXPECT_TRUE(qual.empty());
    EXPECT_EQ(offset, 0);
    EXPECT_TRUE(!ref_offset.has_value());
    EXPECT_TRUE(get<0>(alignment).empty());
    EXPECT_TRUE(get<1>(alignment).empty());
    EXPECT_EQ(flag, 0u);
    EXPECT_EQ(mapq, 0u);
    EXPECT_TRUE(!get<0>(mate).has_value());
    EXPECT_TRUE(!get<1>(mate).has_value());
    EXPECT_EQ(get<2>(mate), int32_t{});
    EXPECT_TRUE(tag_dict.empty());
}

TEST_F(reading_sam, read_in_nothing)
{
    alignment_file_format_sam format;

    std::string file_in_str =
R"(read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read2	42	ref	2	62	1H7M1D1M1S	ref	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1D1M1I1M1I1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)";

    std::istringstream istream(file_in_str);

    alignment_file_header header{};

    for (size_t i = 0; i < 3; ++i)
    {
        ASSERT_NO_THROW(format.read(istream, input_options, std::ignore, header, std::ignore, std::ignore,
                                    std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                    std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore));
    }
}

TEST_F(reading_sam, read_in_alignment_only)
{
    alignment_file_format_sam format;

    std::string file_in_str =
R"(read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read2	42	ref	2	62	1H7M1D1M1S	ref	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1D1M1I1M1I1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)";

    std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>> alignment;
    std::optional<int32_t> ref_id_in;

    std::istringstream istream(file_in_str);

    /*with reference information*/
    for (size_t i = 0; i < 3; ++i)
    {
        ASSERT_NO_THROW(format.read(istream, input_options, ref_sequences, header, std::ignore, std::ignore,
                                    std::ignore, std::ignore, std::ignore, ref_id_in, std::ignore, alignment,
                                    std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore));

        EXPECT_EQ(get<0>(alignment), get<0>(alignments[i]));
        EXPECT_EQ(get<1>(alignment), get<1>(alignments[i]));

        alignment = std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>>{}; // reset
        ref_id_in = 0;

    }

    /*no reference information*/
    using dummy_type = gap_decorator_anchor_set<decltype(ranges::view::repeat_n(dna5{}, size_t{}) |
                                                         std::view::transform(detail::access_restrictor_fn{}))>;
    std::pair<dummy_type, std::vector<gapped<dna5>>> alignment2;
    istream = std::istringstream(file_in_str);

    for (size_t i = 0; i < 3; ++i)
    {
        ASSERT_NO_THROW(format.read(istream, input_options, std::ignore, header, std::ignore, std::ignore, std::ignore,
                                    std::ignore, std::ignore, ref_id_in, std::ignore, alignment2, std::ignore,
                                    std::ignore, std::ignore, std::ignore, std::ignore, std::ignore));

        EXPECT_EQ(get<1>(alignment2), get<1>(alignments[i]));

        alignment2 = std::pair<dummy_type, std::vector<gapped<dna5>>>{}; // reset
        ref_id_in = 0;
    }

    /*no alignment information*/
    std::istringstream istream_empty_cigar(std::string("read1\t41\tref\t1\t61\t*\tref\t10\t300\tACGT\t!##$\n"));

    // with reference sequence information
    ASSERT_NO_THROW(format.read(istream_empty_cigar, input_options, ref_sequences, header, std::ignore, std::ignore,
                                std::ignore, std::ignore, std::ignore, ref_id_in, std::ignore, alignment,
                                std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore));

    EXPECT_TRUE(std::ranges::empty(get<0>(alignment)));
    EXPECT_TRUE(std::ranges::empty(get<1>(alignment)));

    // without reference sequence information
    istream_empty_cigar = std::istringstream(std::string("read1\t41\tref\t1\t61\t*\tref\t10\t300\tACGT\t!##$\n"));
    ASSERT_NO_THROW(format.read(istream_empty_cigar, input_options, std::ignore, header, std::ignore, std::ignore,
                                std::ignore, std::ignore, std::ignore, ref_id_in, std::ignore, alignment2, std::ignore,
                                std::ignore, std::ignore, std::ignore, std::ignore, std::ignore));

    EXPECT_TRUE(std::ranges::empty(get<0>(alignment2)));
    EXPECT_TRUE(std::ranges::empty(get<1>(alignment2)));
}


TEST_F(reading_sam, read_mate_but_not_ref_id)
{
    alignment_file_format_sam format;

    std::string file_in_str =
R"(read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read2	42	ref	2	62	1H7M1D1M1S	=	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1D1M1I1M1I1D1M1S	=	10	300	GGAGTATA	!!*+,-./
)";

    std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> mate;

    std::istringstream istream(file_in_str);

    /*with reference information*/
    for (size_t i = 0; i < 3; ++i)
    {
        ASSERT_NO_THROW(format.read(istream, input_options, ref_sequences, header, std::ignore, std::ignore,
                                    std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                    std::ignore, std::ignore, mate, std::ignore, std::ignore, std::ignore));

        EXPECT_EQ(mate, mates[i]);
        mate = std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>{};
    }

    /*no reference information*/
    istream = std::istringstream(file_in_str);
    for (size_t i = 0; i < 3; ++i)
    {
        ASSERT_NO_THROW(format.read(istream, input_options, std::ignore, header, std::ignore, std::ignore,
                                    std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                    std::ignore, std::ignore, mate, std::ignore, std::ignore, std::ignore));

        EXPECT_EQ(mate, mates[i]);
        mate = std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>{};
    }
}

TEST_F(reading_sam, windows_file)
{
    alignment_file_format_sam format;
    std::string id;
    std::istringstream istream_empty_cigar(std::string("read1\t41\tref\t1\t61\t*\tref\t10\t300\tACGT\t!##$\r\n"));

    // with reference sequence information
    ASSERT_NO_THROW(format.read(istream_empty_cigar, input_options, ref_sequences, header, std::ignore, std::ignore,
                                id, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                std::ignore, std::ignore, std::ignore, std::ignore, std::ignore));

    EXPECT_EQ(id, std::string{"read1"});
}

struct reading_format_error : public reading_sam
{
    std::string seq;
};

TEST_F(reading_format_error, illegal_character_in_seq)
{
    alignment_file_format_sam format;

    std::istringstream istream(std::string("*\t0\t*\t0\t0\t*\t*\t0\t0\tAC!T\t*\n"));

    alignment_file_header header{};
    std::string seq;

    EXPECT_THROW(format.read(istream, input_options, std::ignore, header, seq,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                 format_error);
}

TEST_F(reading_format_error, invalid_arithmetic_value)
{
    alignment_file_format_sam format;

    // invalid value
    std::istringstream istream(std::string("*\t0\t*\t1abc\t0\t*\t*\t0\t0\t*\t*\n"));

    alignment_file_header header{};
    std::optional<int32_t> ref_offset;

    EXPECT_THROW(format.read(istream, input_options, std::ignore, header, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, ref_offset,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                 format_error);

    // overflow error
    istream = std::istringstream(std::string("*\t0\t*\t2147483650\t0\t*\t*\t0\t0\t*\t*\n"));

    EXPECT_THROW(format.read(istream, input_options, std::ignore, header, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, ref_offset,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                 format_error);

    // negative value as ref_offset
    istream = std::istringstream(std::string("*\t0\t*\t-3\t0\t*\t*\t0\t0\t*\t*\n"));

    EXPECT_THROW(format.read(istream, input_options, std::ignore, header, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, ref_offset,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                 format_error);

    // negative value as mate mapping position
    std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> mate;
    istream = std::istringstream(std::string("*\t0\t*\t0\t0\t*\t*\t-3\t0\t*\t*\n"));

    EXPECT_THROW(format.read(istream, input_options, std::ignore, header, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                             mate, std::ignore, std::ignore, std::ignore),
                 format_error);
}

TEST_F(reading_format_error, ref_id_not_in_reference_information)
{
    alignment_file_format_sam format;

    std::istringstream istream(std::string("*\t0\tunknown_ref\t1\t0\t4M\t*\t0\t0\tAAAA\t*\n"));

    std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>> alignment;
    std::optional<int32_t> ref_id_in;

    // with reference information given
    EXPECT_THROW(format.read(istream, input_options,  ref_sequences, header, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, ref_id_in, std::ignore, alignment,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                 format_error);

    istream = std::istringstream(
                  std::string("@HD\tVN:1.6\n@SQ\tSN:ref\tLN:34\n*\t0\tunknown_ref\t1\t0\t4M\t*\t0\t0\tAAAA\t*\n"));

    // with reference information in the header
    using dummy_type = gap_decorator_anchor_set<decltype(ranges::view::repeat_n(dna5{}, size_t{}) |
                                                         std::view::transform(detail::access_restrictor_fn{}))>;
    std::pair<dummy_type, std::vector<gapped<dna5>>> alignment2;
    alignment_file_header<> default_header{};
    EXPECT_THROW(format.read(istream, input_options, std::ignore, default_header, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, ref_id_in, std::ignore, alignment2,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                 format_error);
}

TEST_F(reading_format_error, invalid_sam_tag_format)
{
    alignment_file_format_sam format;

    // type identifier is wrong
    std::istringstream istream(std::string("*\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*\tNM:X:3\n"));

    alignment_file_header header{};
    sam_tag_dictionary dict;

    EXPECT_THROW(format.read(istream, input_options, std::ignore, header, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                             std::ignore, dict, std::ignore, std::ignore),
                 format_error);

    // Array subtype identifier is wrong
    istream = std::istringstream(std::string("*\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*\tNM:B:x3,4\n"));

    EXPECT_THROW(format.read(istream, input_options, std::ignore, header, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                             std::ignore, dict, std::ignore, std::ignore),
                 format_error);
}

TEST_F(reading_format_error, invalid_cigar)
{
    alignment_file_format_sam format;

    // operation P is not supported yet
    std::istringstream istream(std::string("*\t0\t*\t0\t0\t5P\t*\t0\t0\t*\t*\n"));

    alignment_file_header header{};
    std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>> alignment;

    EXPECT_THROW(format.read(istream, input_options, ref_sequences, header, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, alignment, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore),
                 format_error);

    // unkown operation
    istream = std::istringstream(std::string("*\t0\t*\t0\t0\t5Z\t*\t0\t0\t*\t*\n"));
    EXPECT_THROW(format.read(istream, input_options, ref_sequences, header, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, alignment, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore),
                 format_error);

    // negative number as operation count
    istream = std::istringstream(std::string("*\t0\t*\t0\t0\t-5M\t*\t0\t0\t*\t*\n"));
    EXPECT_THROW(format.read(istream, input_options, ref_sequences, header, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, alignment, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore),
                 format_error);

    istream = std::istringstream(std::string("*\t0\t*\t0\t0\t3S4M1I-5M2D2M\t*\t0\t0\t*\t*\n"));
    EXPECT_THROW(format.read(istream, input_options, ref_sequences, header, std::ignore, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore, alignment, std::ignore, std::ignore,
                             std::ignore, std::ignore, std::ignore, std::ignore),
                 format_error);
}

// ----------------------------------------------------------------------------
// writing
// ----------------------------------------------------------------------------

struct writing_sam : public reading_sam
{};

TEST_F(writing_sam, default_options_all_members_specified)
{
    alignment_file_format_sam format;

    std::ostringstream ostream;

    alignment_file_header header{std::vector<std::string>{ref_id}};
    header.ref_id_info.push_back({ref_seq.size(), ""});
    header.ref_dict[ref_id] = 0;

    tag_dicts[0].get<"NM"_tag>() = 7;
    tag_dicts[0].get<"AS"_tag>() = 2;
    tag_dicts[1]["xy"_tag] = std::vector<uint16_t>{3,4,5};

    std::string comp =
R"(@HD	VN:1.6
@SQ	SN:ref	LN:34
read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read2	42	ref	2	62	7M1D1M1S	ref	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1D1M1I1M1I1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)";

    for (size_t i = 0; i < 3; ++i)
        ASSERT_NO_THROW(format.write(ostream, output_options, header, seqs[i], quals[i], ids[i], offsets[i],
                                     std::string{}, 0, ref_offsets[i], alignments[i], flags[i], mapqs[i], mates[i],
                                     tag_dicts[i], 0, 0));

    ostream.flush();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(writing_sam, with_header)
{
    alignment_file_format_sam format;

    std::ostringstream ostream;

    alignment_file_header header{std::vector<std::string>{ref_id}};
    header.sorting = "unknown";
    header.grouping = "none";
    header.ref_id_info.push_back({ref_seq.size(), "AN:other_name"});
    header.ref_dict[ref_id] = 0;
    header.program_infos.push_back({"prog1", "cool_program", "./prog1", "a", "b", "c"});
    header.read_groups.emplace_back("group1", "more info");
    header.comments.push_back("This is a comment.");

    tag_dicts[0].get<"NM"_tag>() = 7;
    tag_dicts[0].get<"AS"_tag>() = 2;
    tag_dicts[1]["xy"_tag] = std::vector<uint16_t>{3,4,5};

    std::string comp =
R"(@HD	VN:1.6	SO:unknown	GO:none
@SQ	SN:ref	LN:34	AN:other_name
@RG	ID:group1	more info
@PG	ID:prog1	PN:cool_program	CL:./prog1	PP:a	DS:b	VN:c
@CO	This is a comment.
read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read2	42	ref	2	62	7M1D1M1S	ref	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1D1M1I1M1I1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)";

    for (size_t i = 0; i < 3; ++i)
        ASSERT_NO_THROW(format.write(ostream, output_options, header, seqs[i], quals[i], ids[i], offsets[i],
                                     std::string{}, 0, ref_offsets[i], alignments[i], flags[i], mapqs[i], mates[i],
                                     tag_dicts[i], 0, 0));

    ostream.flush();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(writing_sam, write_different_header)
{
    std::ostringstream ostream;

    alignment_file_header header{std::vector<std::string>{ref_id}};
    header.ref_id_info.push_back({ref_seq.size(), ""});
    header.ref_dict[ref_id] = 0;

    auto write_header = [&] ()
    {
        alignment_file_format_sam format;

        ASSERT_NO_THROW(format.write(ostream, output_options, header, "", std::vector<phred42>{}, "", 0, "", 0, 0,
                                     std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>>{}, 0, 0,
                                     std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>{},
                                     sam_tag_dictionary{}, 0, 0));
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

TEST_F(writing_sam, special_cases)
{
    alignment_file_format_sam format;

    alignment_file_header header{std::vector<std::string>{ref_id}};
    header.ref_id_info.push_back({ref_seq.size(), ""});
    header.ref_dict[ref_id] = 0;

    std::ostringstream ostream;

    std::string comp =
R"(@HD	VN:1.6
@SQ	SN:ref	LN:34
read1	41	*	1	61	1S1M1D1M1I	*	0	0	ACGT	!##$
)";

    // write an empty std::optional for ref offset and mate
    std::optional<int32_t> rid;
    std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> mate{rid, rid, 0};
    EXPECT_NO_THROW(format.write(ostream, output_options, header, seqs[0], quals[0], ids[0], offsets[0], std::string{},
                                 rid, ref_offsets[0], alignments[0], flags[0],
                                 mapqs[0], mate, tag_dicts[0], 0, 0));
    ostream.flush();
    EXPECT_EQ(ostream.str(), comp);

    ostream = std::ostringstream{};       // clear
    format = alignment_file_format_sam{}; // clear header_was_written

    // write the ref id and mate ref as string
    std::tuple<std::string, std::optional<int32_t>, int32_t> mate_str{"", rid, 0};
    EXPECT_NO_THROW(format.write(ostream, output_options, header, seqs[0], quals[0], ids[0], offsets[0], std::string{},
                                 std::string(""), ref_offsets[0], alignments[0], flags[0],
                                 mapqs[0], mate_str, tag_dicts[0], 0, 0));
    ostream.flush();
    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(writing_sam, format_errors)
{
    alignment_file_format_sam format;

    alignment_file_header header{std::vector<std::string>{ref_id}};
    header.ref_id_info.push_back({ref_seq.size(), ""});
    header.ref_dict[ref_id] = 0;

    std::ostringstream ostream;

    // ensure that only a ref_id that is listed in the header is allowed
    EXPECT_THROW(format.write(ostream, output_options, header, seqs[0], quals[0], ids[0], offsets[0], std::string{},
                              std::string("ref_id_that_does_not_exist"), ref_offsets[0], alignments[0], flags[0],
                              mapqs[0], mates[0], tag_dicts[0], 0, 0),
                 format_error);

    EXPECT_THROW(format.write(ostream, output_options, header, seqs[0], quals[0], ids[0], offsets[0], std::string{},
                              ref_id, -3, alignments[0], flags[0],
                              mapqs[0], mates[0], tag_dicts[0], 0, 0),
                 format_error);
}
