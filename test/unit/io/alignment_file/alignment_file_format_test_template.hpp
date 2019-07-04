// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

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

struct alignment_file_data : public ::testing::Test
{
    alignment_file_data()
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

template <typename format_t>
struct alignment_file_read : public alignment_file_data
{};

TYPED_TEST_CASE_P(alignment_file_read);

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TYPED_TEST_P(alignment_file_read, input_concept)
{
    EXPECT_TRUE((AlignmentFileInputFormat<TypeParam>));
}

// ----------------------------------------------------------------------------
// alignment_file_read
// ----------------------------------------------------------------------------

TYPED_TEST_P(alignment_file_read, header_sucess)
{
    detail::alignment_file_input_format<TypeParam> format;
    typename TestFixture::stream_type istream{this->big_header_input};

    alignment_file_header header{};

    ASSERT_NO_THROW(format.read(istream, input_options, std::ignore, header,  std::ignore, std::ignore, std::ignore,
                                std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                std::ignore, std::ignore, std::ignore, std::ignore, std::ignore));

    EXPECT_EQ(header.format_version, "1.6");
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

TYPED_TEST_P(alignment_file_read, read_in_all_data)
{
    detail::alignment_file_input_format<TypeParam> format;
    typename TestFixture::stream_type istream{this->verbose_reads_input};

    this->tag_dicts[0]["NM"_tag] = -7;
    this->tag_dicts[0]["AS"_tag] = 2;
    this->tag_dicts[0]["CC"_tag] = 300;
    this->tag_dicts[0]["cc"_tag] = -300;
    this->tag_dicts[0]["aa"_tag] = 'c';
    this->tag_dicts[0]["ff"_tag] = 3.1f;
    this->tag_dicts[0]["zz"_tag] = "str";
    this->tag_dicts[1]["bc"_tag] = std::vector<int8_t>{-3};
    this->tag_dicts[1]["bC"_tag] = std::vector<uint8_t>{3u, 200u};
    this->tag_dicts[1]["bs"_tag] = std::vector<int16_t>{-3, 200, -300};
    this->tag_dicts[1]["bS"_tag] = std::vector<uint16_t>{300u, 40u, 500u};
    this->tag_dicts[1]["bi"_tag] = std::vector<int32_t>{-3, 200, -66000};
    this->tag_dicts[1]["bI"_tag] = std::vector<uint32_t>{294967296u};
    this->tag_dicts[1]["bf"_tag] = std::vector<float>{3.5f, 0.1f, 43.8f};

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

    for (size_t i = 0; i < 3; ++i)
    {
        ASSERT_NO_THROW(format.read(istream, input_options, this->ref_sequences, this->header, seq, qual, id, offset,
                                    std::ignore, ref_id_in, ref_offset, alignment, flag, mapq, mate, tag_dict,
                                    std::ignore, std::ignore));

        EXPECT_EQ(seq, this->seqs[i]);
        EXPECT_EQ(id, this->ids[i]);
        EXPECT_EQ(qual, this->quals[i]);
        EXPECT_EQ(offset, this->offsets[i]);
        EXPECT_EQ(ref_id_in, 0);
        EXPECT_EQ(*ref_offset, this->ref_offsets[i]);
        EXPECT_EQ(get<0>(alignment), get<0>(this->alignments[i]));
        EXPECT_EQ(get<1>(alignment), get<1>(this->alignments[i]));
        EXPECT_EQ(flag, this->flags[i]);
        EXPECT_EQ(mapq, this->mapqs[i]);
        EXPECT_EQ(mate, this->mates[i]);
        EXPECT_EQ(tag_dict, this->tag_dicts[i]);

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

TYPED_TEST_P(alignment_file_read, read_in_all_but_empty_data)
{
    detail::alignment_file_input_format<TypeParam> format;
    typename TestFixture::stream_type istream{this->empty_input};

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

    ASSERT_NO_THROW(format.read(istream, input_options, this->ref_sequences, this->header, seq, qual, id, offset, std::ignore,
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

TYPED_TEST_P(alignment_file_read, read_in_nothing)
{
    detail::alignment_file_input_format<TypeParam> format;
    typename TestFixture::stream_type istream{this->simple_three_reads_input};

    alignment_file_header header{};

    for (size_t i = 0; i < 3; ++i)
    {
        ASSERT_NO_THROW(format.read(istream, input_options, std::ignore, header, std::ignore, std::ignore,
                                    std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                    std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore));
    }
}

TYPED_TEST_P(alignment_file_read, read_in_alignment_only_with_ref)
{
    std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>> alignment;
    std::optional<int32_t> ref_id_in;

    {
        detail::alignment_file_input_format<TypeParam> format;
        typename TestFixture::stream_type istream{this->simple_three_reads_input};
        /*with reference information*/
        for (size_t i = 0; i < 3; ++i)
        {
            ASSERT_NO_THROW(format.read(istream, input_options, this->ref_sequences, this->header, std::ignore, std::ignore,
                                        std::ignore, std::ignore, std::ignore, ref_id_in, std::ignore, alignment,
                                        std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore));

            EXPECT_EQ(get<0>(alignment), get<0>(this->alignments[i]));
            EXPECT_EQ(get<1>(alignment), get<1>(this->alignments[i]));

            alignment = std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>>{}; // reset
            ref_id_in = 0;
        }
    }

    {   // empty cigar
        detail::alignment_file_input_format<TypeParam> format;
        typename TestFixture::stream_type istream{this->empty_cigar};
        std::istringstream istream_empty_cigar{};


        ASSERT_NO_THROW(format.read(istream, input_options, this->ref_sequences, this->header, std::ignore, std::ignore,
                                    std::ignore, std::ignore, std::ignore, ref_id_in, std::ignore, alignment,
                                    std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore));

        EXPECT_TRUE(std::ranges::empty(get<0>(alignment)));
        EXPECT_TRUE(std::ranges::empty(get<1>(alignment)));
    }
}

TYPED_TEST_P(alignment_file_read, read_in_alignment_only_without_ref)
{
    using dummy_type = gap_decorator<decltype(view::repeat_n(dna5{}, size_t{}) |
                                                         std::view::transform(detail::access_restrictor_fn{}))>;
    std::optional<int32_t> ref_id_in;
    std::pair<dummy_type, std::vector<gapped<dna5>>> alignment2;

    {
        detail::alignment_file_input_format<TypeParam> format;
        typename TestFixture::stream_type istream{this->simple_three_reads_input};
        alignment_file_header<> default_header{};

        for (size_t i = 0; i < 3; ++i)
        {
            ASSERT_NO_THROW(format.read(istream, input_options, std::ignore, default_header, std::ignore, std::ignore,
                                        std::ignore, std::ignore, std::ignore, ref_id_in, std::ignore, alignment2,
                                        std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore));

            EXPECT_EQ(get<1>(alignment2), get<1>(this->alignments[i]));

            alignment2 = std::pair<dummy_type, std::vector<gapped<dna5>>>{}; // reset
            ref_id_in = 0;
        }
    }

    {   // empty cigar
        detail::alignment_file_input_format<TypeParam> format;
        typename TestFixture::stream_type istream{this->empty_cigar};
        alignment_file_header<> default_header{};

        ASSERT_NO_THROW(format.read(istream, input_options, std::ignore, default_header, std::ignore, std::ignore,
                                    std::ignore, std::ignore, std::ignore, ref_id_in, std::ignore, alignment2,
                                    std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore));

        EXPECT_TRUE(std::ranges::empty(get<0>(alignment2)));
        EXPECT_TRUE(std::ranges::empty(get<1>(alignment2)));
    }
}

TYPED_TEST_P(alignment_file_read, read_mate_but_not_ref_id_with_ref)
{
    std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> mate;

    {   /*with reference information*/
        detail::alignment_file_input_format<TypeParam> format;
        typename TestFixture::stream_type istream{this->simple_three_reads_input};

        for (size_t i = 0; i < 3; ++i)
        {
            ASSERT_NO_THROW(format.read(istream, input_options, this->ref_sequences, this->header, std::ignore,
                                        std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                        std::ignore, std::ignore, std::ignore, mate, std::ignore, std::ignore,
                                        std::ignore));

            EXPECT_EQ(mate, this->mates[i]);
            mate = std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>{};
        }
    }
}

TYPED_TEST_P(alignment_file_read, read_mate_but_not_ref_id_without_ref)
{
    std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> mate;

    {   /*no reference information*/
        detail::alignment_file_input_format<TypeParam> format;
        typename TestFixture::stream_type istream{this->simple_three_reads_input};
        alignment_file_header<> default_header{};

        for (size_t i = 0; i < 3; ++i)
        {
            ASSERT_NO_THROW(format.read(istream, input_options, std::ignore, default_header, std::ignore, std::ignore,
                                        std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore,
                                        std::ignore, std::ignore, mate, std::ignore, std::ignore, std::ignore));

            EXPECT_EQ(mate, this->mates[i]);
            mate = std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>{};
        }
    }
}

TYPED_TEST_P(alignment_file_read, format_error_ref_id_not_in_reference_information)
{
    using dummy_type = gap_decorator<decltype(view::repeat_n(dna5{}, size_t{}) |
                                                         std::view::transform(detail::access_restrictor_fn{}))>;
    std::optional<int32_t> ref_id_in;
    std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>> alignment;
    std::pair<dummy_type, std::vector<gapped<dna5>>> alignment2;

    {   // with reference information given
        detail::alignment_file_input_format<TypeParam> format;
        typename TestFixture::stream_type istream{this->unknown_ref};

        EXPECT_THROW(format.read(istream, input_options,  this->ref_sequences, this->header, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore, ref_id_in, std::ignore, alignment,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                     format_error);
    }

    {   // with reference information in the header
        detail::alignment_file_input_format<TypeParam> format;
        typename TestFixture::stream_type istream{this->unknown_ref_header};
        alignment_file_header<> default_header{};

        EXPECT_THROW(format.read(istream, input_options, std::ignore, default_header, std::ignore, std::ignore,
                                 std::ignore, std::ignore, std::ignore, ref_id_in, std::ignore, alignment2,
                                 std::ignore, std::ignore, std::ignore, std::ignore, std::ignore, std::ignore),
                     format_error);
    }
}

// ----------------------------------------------------------------------------
// alignment_file_write
// ----------------------------------------------------------------------------

template <typename format_type>
struct alignment_file_write : public alignment_file_read<format_type>
{};

TYPED_TEST_CASE_P(alignment_file_write);

TYPED_TEST_P(alignment_file_write, output_concept)
{
    EXPECT_TRUE((AlignmentFileOutputFormat<TypeParam>));
}

TYPED_TEST_P(alignment_file_write, write_empty_members)
{
    detail::alignment_file_output_format<TypeParam> format;

    std::ostringstream ostream;
    {
        alignment_file_header header{std::vector<std::string>{this->ref_id}};
        header.ref_id_info.push_back({this->ref_seq.size(), ""});
        header.ref_dict[this->ref_id] = 0;

        using default_align_t = std::pair<std::span<gapped<char>>, std::span<gapped<char>>>;
        using default_mate_t  = std::tuple<std::string_view, std::optional<int32_t>, int32_t>;

        ASSERT_NO_THROW(format.write(ostream, output_options, header, std::string_view{}, std::string_view{},
                                     std::string_view{}, 0, std::string_view{}, std::string_view{},
                                     std::optional<int32_t>{std::nullopt}, default_align_t{}, 0, 0,
                                     default_mate_t{}, sam_tag_dictionary{}, 0, 0));
    }

    ostream.flush();
    EXPECT_EQ(ostream.str(), this->empty_input);
}

TYPED_TEST_P(alignment_file_write, default_options_all_members_specified)
{
    detail::alignment_file_output_format<TypeParam> format;

    std::ostringstream ostream;

    alignment_file_header header{std::vector<std::string>{this->ref_id}};
    header.ref_id_info.push_back({this->ref_seq.size(), ""});
    header.ref_dict[this->ref_id] = 0;

    this->tag_dicts[0]["NM"_tag] = 7;
    this->tag_dicts[0]["AS"_tag] = 2;
    this->tag_dicts[1]["xy"_tag] = std::vector<uint16_t>{3,4,5};

    for (size_t i = 0; i < 3; ++i)
        ASSERT_NO_THROW(format.write(ostream, output_options, header, this->seqs[i], this->quals[i], this->ids[i],
                                     this->offsets[i], std::string{}, 0, this->ref_offsets[i], this->alignments[i],
                                     this->flags[i], this->mapqs[i], this->mates[i], this->tag_dicts[i], 0, 0));

    ostream.flush();

    EXPECT_EQ(ostream.str(), this->simple_three_reads_output);
}

TYPED_TEST_P(alignment_file_write, write_ref_id_with_different_types)
{
    detail::alignment_file_output_format<TypeParam> format;

    std::ostringstream ostream;

    alignment_file_header header{std::vector<std::string>{this->ref_id}};
    header.ref_id_info.push_back({this->ref_seq.size(), ""});
    header.ref_dict[this->ref_id] = 0;

    this->tag_dicts[0]["NM"_tag] = 7;
    this->tag_dicts[0]["AS"_tag] = 2;
    this->tag_dicts[1]["xy"_tag] = std::vector<uint16_t>{3,4,5};

    // header ref_id_type is std::string

    // std::string
    ASSERT_NO_THROW(format.write(ostream, output_options, header, this->seqs[0], this->quals[0], this->ids[0],
                                 this->offsets[0], std::string{},
    /*----------------------->*/ this->ref_id,
                                 this->ref_offsets[0], this->alignments[0],
                                 this->flags[0], this->mapqs[0], this->mates[0], this->tag_dicts[0], 0, 0));
    // std::string_view
    ASSERT_NO_THROW(format.write(ostream, output_options, header, this->seqs[1], this->quals[1], this->ids[1],
                                 this->offsets[1], std::string{},
    /*----------------------->*/ std::string_view{this->ref_id},
                                 this->ref_offsets[1], this->alignments[1],
                                 this->flags[1], this->mapqs[1], this->mates[1], this->tag_dicts[1], 0, 0));

    // view on string
    ASSERT_NO_THROW(format.write(ostream, output_options, header, this->seqs[2], this->quals[2], this->ids[2],
                                 this->offsets[2], std::string{},
    /*----------------------->*/ this->ref_id | view::take(20),
                                 this->ref_offsets[2], this->alignments[2],
                                 this->flags[2], this->mapqs[2], this->mates[2], this->tag_dicts[2], 0, 0));

    ostream.flush();

    EXPECT_EQ(ostream.str(), this->simple_three_reads_output);
}

TYPED_TEST_P(alignment_file_write, with_header)
{
    detail::alignment_file_output_format<TypeParam> format;

    std::ostringstream ostream;

    alignment_file_header header{std::vector<std::string>{this->ref_id}};
    header.sorting = "unknown";
    header.grouping = "none";
    header.ref_id_info.push_back({this->ref_seq.size(), "AN:other_name"});
    header.ref_dict[this->ref_id] = 0;
    header.program_infos.push_back({"prog1", "cool_program", "./prog1", "a", "b", "c"});
    header.read_groups.emplace_back("group1", "more info");
    header.comments.push_back("This is a comment.");

    this->tag_dicts[0]["NM"_tag] = -7;
    this->tag_dicts[0]["AS"_tag] = 2;
    this->tag_dicts[0]["CC"_tag] = 300;
    this->tag_dicts[0]["cc"_tag] = -300;
    this->tag_dicts[0]["aa"_tag] = 'c';
    this->tag_dicts[0]["ff"_tag] = 3.1f;
    this->tag_dicts[0]["zz"_tag] = "str";
    this->tag_dicts[1]["bc"_tag] = std::vector<int8_t>{-3};
    this->tag_dicts[1]["bC"_tag] = std::vector<uint8_t>{3u, 200u};
    this->tag_dicts[1]["bs"_tag] = std::vector<int16_t>{-3, 200, -300};
    this->tag_dicts[1]["bS"_tag] = std::vector<uint16_t>{300u, 40u, 500u};
    this->tag_dicts[1]["bi"_tag] = std::vector<int32_t>{-3, 200, -66000};
    this->tag_dicts[1]["bI"_tag] = std::vector<uint32_t>{294967296u};
    this->tag_dicts[1]["bf"_tag] = std::vector<float>{3.5f, 0.1f, 43.8f};

    for (size_t i = 0; i < 3; ++i)
        ASSERT_NO_THROW(format.write(ostream, output_options, header, this->seqs[i], this->quals[i], this->ids[i],
        							 this->offsets[i], std::string{}, 0, this->ref_offsets[i], this->alignments[i],
        							 this->flags[i], this->mapqs[i], this->mates[i], this->tag_dicts[i], 0, 0));

    ostream.flush();

    EXPECT_EQ(ostream.str(), this->verbose_output);
}

TYPED_TEST_P(alignment_file_write, special_cases)
{
    detail::alignment_file_output_format<TypeParam> format;

    alignment_file_header header{std::vector<std::string>{this->ref_id}};
    header.ref_id_info.push_back({this->ref_seq.size(), ""});
    header.ref_dict[this->ref_id] = 0;

    std::ostringstream ostream;

    // write an empty std::optional for ref offset and mate
    std::optional<int32_t> rid;
    std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> mate{rid, rid, 0};

    EXPECT_NO_THROW(format.write(ostream, output_options, header, this->seqs[0], this->quals[0], this->ids[0],
                                 this->offsets[0], std::string{}, rid, this->ref_offsets[0], this->alignments[0],
                                 this->flags[0], this->mapqs[0], mate, this->tag_dicts[0], 0, 0));
    ostream.flush();
    EXPECT_EQ(ostream.str(), this->special_output);

    ostream = std::ostringstream{}; // clear
    format = detail::alignment_file_output_format<TypeParam>{}; // clear header_was_written

    // write the ref id and mate ref as string
    std::tuple<std::string, std::optional<int32_t>, int32_t> mate_str{"", rid, 0};

    /*EXPECT_NO_THROW(*/format.write(ostream, output_options, header, this->seqs[0], this->quals[0], this->ids[0],
                                 this->offsets[0], std::string{}, std::string(""), this->ref_offsets[0],
                                 this->alignments[0], this->flags[0], this->mapqs[0], mate_str,
                                 this->tag_dicts[0], 0, 0)/*)*/;
    ostream.flush();
    EXPECT_EQ(ostream.str(), this->special_output);
}

TYPED_TEST_P(alignment_file_write, format_errors)
{
    detail::alignment_file_output_format<TypeParam> format;

    alignment_file_header header{std::vector<std::string>{this->ref_id}};
    header.ref_id_info.push_back({this->ref_seq.size(), ""});
    header.ref_dict[this->ref_id] = 0;

    std::ostringstream ostream;

    // ensure that only a ref_id that is listed in the header is allowed
    EXPECT_THROW(format.write(ostream, output_options, header, this->seqs[0], this->quals[0], this->ids[0],
                              this->offsets[0], std::string{}, std::string("ref_id_that_does_not_exist"),
                              this->ref_offsets[0], this->alignments[0], this->flags[0],
                              this->mapqs[0], this->mates[0], this->tag_dicts[0], 0, 0),
                 format_error);

    EXPECT_THROW(format.write(ostream, output_options, header, this->seqs[0], this->quals[0], this->ids[0],
    						  this->offsets[0], std::string{}, this->ref_id, -3, this->alignments[0], this->flags[0],
                              this->mapqs[0], this->mates[0], this->tag_dicts[0], 0, 0),
                 format_error);
}

REGISTER_TYPED_TEST_CASE_P(alignment_file_read,
                           input_concept,
                           header_sucess,
                           read_in_all_data,
                           read_in_all_but_empty_data,
                           read_in_nothing,
                           read_in_alignment_only_with_ref,
                           read_in_alignment_only_without_ref,
                           read_mate_but_not_ref_id_with_ref,
                           read_mate_but_not_ref_id_without_ref,
                           format_error_ref_id_not_in_reference_information);

REGISTER_TYPED_TEST_CASE_P(alignment_file_write,
                           write_empty_members,
                           output_concept,
                           default_options_all_members_specified,
                           write_ref_id_with_different_types,
                           with_header,
                           special_cases,
                           format_errors);
