// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/alignment_file/input_format_concept.hpp>
#include <seqan3/io/alignment_file/output_format_concept.hpp>
#include <seqan3/io/alignment_file/input.hpp>
#include <seqan3/io/alignment_file/output.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/pretty_printing.hpp>

using seqan3::operator""_cigar_op;
using seqan3::operator""_dna5;
using seqan3::operator""_phred42;
using seqan3::operator""_tag;

// global variables for reuse
seqan3::alignment_file_input_options<seqan3::dna5> input_options;
seqan3::alignment_file_output_options output_options;

struct alignment_file_data : public ::testing::Test
{
    alignment_file_data()
    {
        ref_sequences = std::vector<seqan3::dna5_vector>{ref_seq};
        ref_ids = std::vector<std::string>{ref_id};
        header = seqan3::alignment_file_header{ref_ids};
        header.ref_id_info.emplace_back(ref_seq.size(), "");
        header.ref_dict[header.ref_ids()[0]] = 0; // set up header which is otherwise done on file level
    }

    std::vector<seqan3::dna5_vector> seqs
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

    std::vector<std::vector<seqan3::phred42>> quals
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

    seqan3::dna5_vector ref_seq = "ACTGATCGAGAGGATCTAGAGGAGATCGTAGGAC"_dna5;

    std::vector<seqan3::gapped<seqan3::dna5>> ref_seq_gapped1 = {'A'_dna5, 'C'_dna5, 'T'_dna5, seqan3::gap{}};
    std::vector<seqan3::gapped<seqan3::dna5>> ref_seq_gapped2 = {'C'_dna5, 'T'_dna5, 'G'_dna5, 'A'_dna5,
                                                                 'T'_dna5, 'C'_dna5, 'G'_dna5, 'A'_dna5, 'G'_dna5};
    std::vector<seqan3::gapped<seqan3::dna5>> ref_seq_gapped3 = {'T'_dna5, seqan3::gap{}, 'G'_dna5, seqan3::gap{},
                                                                 'A'_dna5, seqan3::gap{}, 'T'_dna5, 'C'_dna5};

    std::string ref_id = "ref";

    std::vector<int32_t> ref_offsets
    {
        0,
        1,
        2
    };

    std::vector<std::pair<std::vector<seqan3::gapped<seqan3::dna5>>,
                          std::vector<seqan3::gapped<seqan3::dna5>>>> alignments
    {
        {ref_seq_gapped1, std::vector<seqan3::gapped<seqan3::dna5>>{'C'_dna5, seqan3::gap{}, 'G'_dna5, 'T'_dna5}},
        {ref_seq_gapped2, std::vector<seqan3::gapped<seqan3::dna5>>{'A'_dna5, 'G'_dna5, 'G'_dna5, 'C'_dna5, 'T'_dna5,
                                                                    'G'_dna5, 'N'_dna5, seqan3::gap{}, 'A'_dna5}},
        {ref_seq_gapped3, std::vector<seqan3::gapped<seqan3::dna5>>{'G'_dna5, seqan3::gap{}, 'A'_dna5, 'G'_dna5,
                                                                    'T'_dna5, 'A'_dna5, seqan3::gap{}, 'T'_dna5}}
    };

    std::vector<seqan3::sam_flag> flags
    {
        seqan3::sam_flag{41u},
        seqan3::sam_flag{42u},
        seqan3::sam_flag{43u}
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

    std::vector<seqan3::sam_tag_dictionary> tag_dicts
    {
        seqan3::sam_tag_dictionary{},
        seqan3::sam_tag_dictionary{},
        seqan3::sam_tag_dictionary{}
    };

    std::vector<seqan3::dna5_vector> ref_sequences{};
    std::vector<std::string> ref_ids{};
    seqan3::alignment_file_header<std::vector<std::string>> header{};
};

template <typename format_t>
struct alignment_file_read : public alignment_file_data
{};

TYPED_TEST_SUITE_P(alignment_file_read);

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TYPED_TEST_P(alignment_file_read, input_concept)
{
    EXPECT_TRUE((seqan3::alignment_file_input_format<TypeParam>));
}

// ----------------------------------------------------------------------------
// alignment_file_read
// ----------------------------------------------------------------------------

TYPED_TEST_P(alignment_file_read, header_sucess)
{
    typename TestFixture::stream_type istream{this->big_header_input};

    seqan3::alignment_file_input fin{istream, TypeParam{}};
    auto & header = fin.header();

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
    typename TestFixture::stream_type istream{this->verbose_reads_input};
    seqan3::alignment_file_input fin{istream, this->ref_ids, this->ref_sequences, TypeParam{}};

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

    size_t i{0};
    for (auto & rec : fin)
    {
        EXPECT_EQ(seqan3::get<seqan3::field::seq>(rec), this->seqs[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::id>(rec), this->ids[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::qual>(rec), this->quals[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::offset>(rec), this->offsets[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::ref_id>(rec), 0);
        EXPECT_EQ(*seqan3::get<seqan3::field::ref_offset>(rec), this->ref_offsets[i]);
        EXPECT_RANGE_EQ(std::get<0>(seqan3::get<seqan3::field::alignment>(rec)), std::get<0>(this->alignments[i]));
        EXPECT_RANGE_EQ(std::get<1>(seqan3::get<seqan3::field::alignment>(rec)), std::get<1>(this->alignments[i]));
        EXPECT_EQ(seqan3::get<seqan3::field::flag>(rec), this->flags[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::mapq>(rec), this->mapqs[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::mate>(rec), this->mates[i]);
        EXPECT_EQ(seqan3::get<seqan3::field::tags>(rec), this->tag_dicts[i]);
        ++i;
    }
}

TYPED_TEST_P(alignment_file_read, read_in_all_but_empty_data)
{
    typename TestFixture::stream_type istream{this->empty_input};
    seqan3::alignment_file_input fin{istream, this->ref_ids, this->ref_sequences, TypeParam{}};

    EXPECT_TRUE(seqan3::get<seqan3::field::seq>(*fin.begin()).empty());
    EXPECT_TRUE(seqan3::get<seqan3::field::id>(*fin.begin()).empty());
    EXPECT_TRUE(seqan3::get<seqan3::field::qual>(*fin.begin()).empty());
    EXPECT_EQ(seqan3::get<seqan3::field::offset>(*fin.begin()), 0);
    EXPECT_TRUE(!seqan3::get<seqan3::field::ref_id>(*fin.begin()).has_value());
    EXPECT_TRUE(!seqan3::get<seqan3::field::ref_offset>(*fin.begin()).has_value());
    EXPECT_TRUE(std::ranges::empty(std::get<0>(seqan3::get<seqan3::field::alignment>(*fin.begin()))));
    EXPECT_TRUE(std::ranges::empty(std::get<1>(seqan3::get<seqan3::field::alignment>(*fin.begin()))));
    EXPECT_EQ(seqan3::get<seqan3::field::flag>(*fin.begin()), seqan3::sam_flag{0u});
    EXPECT_EQ(seqan3::get<seqan3::field::mapq>(*fin.begin()), 0u);
    EXPECT_TRUE(!std::get<0>(seqan3::get<seqan3::field::mate>(*fin.begin())).has_value());
    EXPECT_TRUE(!std::get<1>(seqan3::get<seqan3::field::mate>(*fin.begin())).has_value());
    EXPECT_EQ(std::get<2>(seqan3::get<seqan3::field::mate>(*fin.begin())), int32_t{});
    EXPECT_TRUE(seqan3::get<seqan3::field::tags>(*fin.begin()).empty());
}

TYPED_TEST_P(alignment_file_read, read_in_almost_nothing)
{
    typename TestFixture::stream_type istream{this->simple_three_reads_input};
    seqan3::alignment_file_input fin{istream, TypeParam{}, seqan3::fields<seqan3::field::mapq>{}};

    size_t i{0};
    for (auto & [mapq] : fin)
        EXPECT_EQ(mapq, this->mapqs[i++]);
}

TYPED_TEST_P(alignment_file_read, read_in_alignment_only_with_ref)
{
    {
        typename TestFixture::stream_type istream{this->simple_three_reads_input};
        seqan3::alignment_file_input fin{istream,
                                         this->ref_ids,
                                         this->ref_sequences,
                                         TypeParam{},
                                         seqan3::fields<seqan3::field::alignment>{}};

        size_t i{0};
        for (auto & [alignment] : fin)
        {
            EXPECT_RANGE_EQ(std::get<0>(alignment), std::get<0>(this->alignments[i]));
            EXPECT_RANGE_EQ(std::get<1>(alignment), std::get<1>(this->alignments[i]));
            ++i;
        }
    }

    {   // empty cigar
        typename TestFixture::stream_type istream{this->empty_cigar};
        seqan3::alignment_file_input fin{istream,
                                         this->ref_ids,
                                         this->ref_sequences,
                                         TypeParam{},
                                         seqan3::fields<seqan3::field::alignment>{}};

        EXPECT_TRUE(std::ranges::empty(std::get<0>(seqan3::get<seqan3::field::alignment>(*fin.begin()))));
        EXPECT_TRUE(std::ranges::empty(std::get<1>(seqan3::get<seqan3::field::alignment>(*fin.begin()))));
    }
}

TYPED_TEST_P(alignment_file_read, read_in_alignment_only_without_ref)
{
    {
        typename TestFixture::stream_type istream{this->simple_three_reads_input};
        seqan3::alignment_file_input fin{istream, TypeParam{}, seqan3::fields<seqan3::field::alignment>{}};

        size_t i{0};
        for (auto & [alignment] : fin)
        {
            EXPECT_RANGE_EQ(std::get<1>(alignment), std::get<1>(this->alignments[i++]));
            auto & ref_aln = std::get<0>(alignment);
            EXPECT_THROW((ref_aln[0]), std::logic_error); // access on a dummy seq is not allowed
        }
    }

    {   // empty cigar
        typename TestFixture::stream_type istream{this->empty_cigar};
        seqan3::alignment_file_input fin{istream, TypeParam{}, seqan3::fields<seqan3::field::alignment>{}};

        EXPECT_TRUE(std::ranges::empty(std::get<0>(seqan3::get<seqan3::field::alignment>(*fin.begin()))));
        EXPECT_TRUE(std::ranges::empty(std::get<1>(seqan3::get<seqan3::field::alignment>(*fin.begin()))));
    }
}

TYPED_TEST_P(alignment_file_read, read_mate_but_not_ref_id_with_ref)
{
    {   /*with reference information*/
        typename TestFixture::stream_type istream{this->simple_three_reads_input};
        seqan3::alignment_file_input fin{istream,
                                         this->ref_ids,
                                         this->ref_sequences,
                                         TypeParam{},
                                         seqan3::fields<seqan3::field::mate>{}};

        size_t i{0};
        for (auto & [mate] : fin)
            EXPECT_EQ(mate, this->mates[i++]);
    }
}

TYPED_TEST_P(alignment_file_read, read_mate_but_not_ref_id_without_ref)
{
    std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> mate;

    {   /*no reference information*/
        typename TestFixture::stream_type istream{this->simple_three_reads_input};
        seqan3::alignment_file_input fin{istream, TypeParam{}, seqan3::fields<seqan3::field::mate>{}};

        size_t i{0};
        for (auto & [mate] : fin)
            EXPECT_EQ(mate, this->mates[i++]);
    }
}

TYPED_TEST_P(alignment_file_read, cigar_vector)
{
    std::vector<std::vector<seqan3::cigar>> expected
    {
        {{1, 'S'_cigar_op}, {1, 'M'_cigar_op}, {1, 'D'_cigar_op}, {1, 'M'_cigar_op}, {1, 'I'_cigar_op}},
        {{1, 'H'_cigar_op}, {7, 'M'_cigar_op}, {1, 'D'_cigar_op}, {1, 'M'_cigar_op}, {1, 'S'_cigar_op},
         {2, 'H'_cigar_op}},
        {{1, 'S'_cigar_op}, {1, 'M'_cigar_op}, {1, 'P'_cigar_op}, {1, 'M'_cigar_op}, {1, 'I'_cigar_op},
         {1, 'M'_cigar_op}, {1, 'I'_cigar_op}, {1, 'D'_cigar_op}, {1, 'M'_cigar_op}, {1, 'S'_cigar_op}}
    };

    typename TestFixture::stream_type istream{this->simple_three_reads_input};
    seqan3::alignment_file_input fin{istream, TypeParam{}, seqan3::fields<seqan3::field::cigar>{}};

    size_t i{0};
    for (auto & [cigar_v] : fin)
        EXPECT_EQ(cigar_v, expected[i++]);
}

TYPED_TEST_P(alignment_file_read, format_error_ref_id_not_in_reference_information)
{
    {   // with reference information given
        typename TestFixture::stream_type istream{this->unknown_ref};
        seqan3::alignment_file_input fin{istream, this->ref_ids, this->ref_sequences, TypeParam{}};
        EXPECT_THROW((fin.begin()), seqan3::format_error);
    }

    {   // with reference information in the header
        typename TestFixture::stream_type istream{this->unknown_ref_header};
        seqan3::alignment_file_input fin{istream, TypeParam{}};
        EXPECT_THROW((fin.begin()), seqan3::format_error);
    }
}

// ----------------------------------------------------------------------------
// alignment_file_write
// ----------------------------------------------------------------------------

// Note that these differ from the alignment_file_output default fields:
// 1. They don't contain field::bit_score and field::evalue since these belong to the BLAST format.
// 2. field::alignment and field::cigar are redundant. Since field::alignment is the more complex one it is chosen here.
//    The behaviour if both are given is tested in a separate test.
using sam_fields = seqan3::fields<seqan3::field::header_ptr,
                                  seqan3::field::id,
                                  seqan3::field::flag,
                                  seqan3::field::ref_id,
                                  seqan3::field::ref_offset,
                                  seqan3::field::mapq,
                                  seqan3::field::alignment,
                                  seqan3::field::offset,
                                  seqan3::field::mate,
                                  seqan3::field::seq,
                                  seqan3::field::qual,
                                  seqan3::field::tags>;

template <typename format_type>
struct alignment_file_write : public alignment_file_read<format_type>
{
    std::ostringstream ostream;
};

TYPED_TEST_SUITE_P(alignment_file_write);

TYPED_TEST_P(alignment_file_write, output_concept)
{
    EXPECT_TRUE((seqan3::alignment_file_output_format<TypeParam>));
}

TYPED_TEST_P(alignment_file_write, write_empty_members)
{
    {
        seqan3::alignment_file_output fout{this->ostream, TypeParam{}, sam_fields{}};

        using default_align_t = std::pair<std::span<seqan3::gapped<char>>, std::span<seqan3::gapped<char>>>;
        using default_mate_t  = std::tuple<std::string_view, std::optional<int32_t>, int32_t>;

        fout.emplace_back(&(this->header),
                          std::string_view{},
                          seqan3::sam_flag::none,
                          std::string_view{},
                          -1,
                          0,
                          default_align_t{},
                          0,
                          default_mate_t{},
                          std::string_view{},
                          std::string_view{},
                          seqan3::sam_tag_dictionary{});
    }

    this->ostream.flush();
    EXPECT_EQ(this->ostream.str(), this->empty_input);
}

TYPED_TEST_P(alignment_file_write, default_options_all_members_specified)
{
    this->tag_dicts[0]["NM"_tag] = 7;
    this->tag_dicts[0]["AS"_tag] = 2;
    this->tag_dicts[1]["xy"_tag] = std::vector<uint16_t>{3,4,5};

    {
        seqan3::alignment_file_output fout{this->ostream, TypeParam{}, sam_fields{}};

        for (size_t i = 0; i < 3; ++i)
        {
            ASSERT_NO_THROW(fout.emplace_back(&(this->header), this->ids[i], this->flags[i], 0/*ref_id*/,
                                              this->ref_offsets[i], this->mapqs[i], this->alignments[i],
                                              this->offsets[i], this->mates[i], this->seqs[i], this->quals[i],
                                              this->tag_dicts[i]));
        }
    }
    this->ostream.flush();

    EXPECT_EQ(this->ostream.str(), this->simple_three_reads_output);
}

TYPED_TEST_P(alignment_file_write, write_ref_id_with_different_types)
{
    this->tag_dicts[0]["NM"_tag] = 7;
    this->tag_dicts[0]["AS"_tag] = 2;
    this->tag_dicts[1]["xy"_tag] = std::vector<uint16_t>{3,4,5};

    {
        // header ref_id_type is std::string
        seqan3::alignment_file_output fout{this->ostream, TypeParam{}, sam_fields{}};

        // std::string
        ASSERT_NO_THROW(fout.emplace_back(&(this->header), this->ids[0], this->flags[0],
        /*----------------------->*/      this->ref_id,
                                          this->ref_offsets[0], this->mapqs[0], this->alignments[0], this->offsets[0],
                                          this->mates[0], this->seqs[0], this->quals[0], this->tag_dicts[0]));

        // std::string_view
        ASSERT_NO_THROW(fout.emplace_back(&(this->header), this->ids[1], this->flags[1],
        /*----------------------->*/      std::string_view{this->ref_id},
                                          this->ref_offsets[1], this->mapqs[1], this->alignments[1], this->offsets[1],
                                          this->mates[1], this->seqs[1], this->quals[1], this->tag_dicts[1]));

        // view on string
        ASSERT_NO_THROW(fout.emplace_back(&(this->header), this->ids[2], this->flags[2],
        /*----------------------->*/      this->ref_id | seqan3::views::take(20),
                                          this->ref_offsets[2], this->mapqs[2], this->alignments[2], this->offsets[2],
                                          this->mates[2], this->seqs[2], this->quals[2], this->tag_dicts[2]));
    }

    this->ostream.flush();

    EXPECT_EQ(this->ostream.str(), this->simple_three_reads_output);
}

TYPED_TEST_P(alignment_file_write, with_header)
{
    seqan3::alignment_file_header header{std::vector<std::string>{this->ref_id}};
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

    {
        seqan3::alignment_file_output fout{this->ostream, TypeParam{}, sam_fields{}};

        for (size_t i = 0; i < 3; ++i)
        {
            ASSERT_NO_THROW(fout.emplace_back(&header, this->ids[i], this->flags[i], 0/*ref_id*/,
                                              this->ref_offsets[i], this->mapqs[i], this->alignments[i],
                                              this->offsets[i], this->mates[i], this->seqs[i], this->quals[i],
                                              this->tag_dicts[i]));
        }
    }

    this->ostream.flush();

    EXPECT_EQ(this->ostream.str(), this->verbose_output);
}

TYPED_TEST_P(alignment_file_write, cigar_vector)
{
    std::vector<std::vector<seqan3::cigar>> cigar_v
    {
        {{1, 'S'_cigar_op}, {1, 'M'_cigar_op}, {1, 'D'_cigar_op}, {1, 'M'_cigar_op}, {1, 'I'_cigar_op}},
        {{1, 'H'_cigar_op}, {7, 'M'_cigar_op}, {1, 'D'_cigar_op}, {1, 'M'_cigar_op}, {1, 'S'_cigar_op},
         {2, 'H'_cigar_op}},
        {{1, 'S'_cigar_op}, {1, 'M'_cigar_op}, {1, 'P'_cigar_op}, {1, 'M'_cigar_op}, {1, 'I'_cigar_op},
         {1, 'M'_cigar_op}, {1, 'I'_cigar_op}, {1, 'D'_cigar_op}, {1, 'M'_cigar_op}, {1, 'S'_cigar_op}}
    };

    this->tag_dicts[0]["NM"_tag] = 7;
    this->tag_dicts[0]["AS"_tag] = 2;
    this->tag_dicts[1]["xy"_tag] = std::vector<uint16_t>{3,4,5};

    {
        seqan3::alignment_file_output fout{this->ostream, TypeParam{}}; // default fields contain CIGAR and alignment
        for (size_t i = 0ul; i < 3ul; ++i)
        {
            ASSERT_NO_THROW(fout.emplace_back(this->seqs[i],
                                              this->ids[i],
                                              this->offsets[i],
                                              this->ref_seq,
                                              0/*ref_id*/,
                                              this->ref_offsets[i],
                                              this->alignments[i],
                                              cigar_v[i],
                                              this->mapqs[i],
                                              this->quals[i],
                                              this->flags[i],
                                              this->mates[i],
                                              this->tag_dicts[i],
                                              0/*evalue*/,
                                              0/*bitscore*/,
                                              &(this->header)));
        }
    }

    this->ostream.flush();
    // compare to original input because hard clipping is preserved when writing the cigar vector directly
    EXPECT_EQ(this->ostream.str(), this->simple_three_reads_input);

    this->ostream = std::ostringstream{}; // clear

    // 2. Write only the cigar, not the alignment
    {
        seqan3::alignment_file_output fout{this->ostream, TypeParam{}, seqan3::fields<seqan3::field::header_ptr,
                                                                                      seqan3::field::id,
                                                                                      seqan3::field::flag,
                                                                                      seqan3::field::ref_id,
                                                                                      seqan3::field::ref_offset,
                                                                                      seqan3::field::mapq,
                                                                                      seqan3::field::cigar,
                                                                                           // cigar instead of alignment
                                                                                      seqan3::field::offset,
                                                                                      seqan3::field::mate,
                                                                                      seqan3::field::seq,
                                                                                      seqan3::field::qual,
                                                                                      seqan3::field::tags>{}};

        for (size_t i = 0ul; i < 3ul; ++i)
        {
            ASSERT_NO_THROW(fout.emplace_back(&(this->header), this->ids[i], this->flags[i], 0/*ref_id*/,
                                              this->ref_offsets[i], this->mapqs[i], cigar_v[i],
                                              this->offsets[i], this->mates[i], this->seqs[i], this->quals[i],
                                              this->tag_dicts[i]));
        }
    }

    this->ostream.flush();
    // compare to original input because hard clipping is preserved when writing the cigar vector directly
    EXPECT_EQ(this->ostream.str(), this->simple_three_reads_input);
}

TYPED_TEST_P(alignment_file_write, special_cases)
{
    std::optional<int32_t> rid;

    // write an empty std::optional for ref id and mate
    {
        std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> mate{rid, rid, 0};

        seqan3::alignment_file_output fout{this->ostream, TypeParam{}, sam_fields{}};

        // std::string
        ASSERT_NO_THROW(fout.emplace_back(&(this->header), this->ids[0], this->flags[0], rid,
                                          this->ref_offsets[0], this->mapqs[0], this->alignments[0], this->offsets[0],
                                          mate, this->seqs[0], this->quals[0], this->tag_dicts[0]));

    }

    this->ostream.flush();
    EXPECT_EQ(this->ostream.str(), this->special_output);

    this->ostream = std::ostringstream{}; // clear

    {
        // write the ref id and mate ref as string
        std::tuple<std::string, std::optional<int32_t>, int32_t> mate_str{"", rid, 0};

        seqan3::alignment_file_output fout{this->ostream, TypeParam{}, sam_fields{}};

        // std::string
        ASSERT_NO_THROW(fout.emplace_back(&(this->header), this->ids[0], this->flags[0], std::string(""),
                                          this->ref_offsets[0], this->mapqs[0], this->alignments[0], this->offsets[0],
                                          mate_str, this->seqs[0], this->quals[0], this->tag_dicts[0]));

    }

    this->ostream.flush();
    EXPECT_EQ(this->ostream.str(), this->special_output);
}

TYPED_TEST_P(alignment_file_write, format_errors)
{
    seqan3::alignment_file_output fout{this->ostream, TypeParam{}, sam_fields{}};

    // ensure that only a ref_id that is listed in the header is allowed
    EXPECT_THROW(fout.emplace_back(&(this->header), this->ids[0], this->flags[0],
                                   std::string("ref_id_that_does_not_exist"),
                                   this->ref_offsets[0], this->mapqs[0], this->alignments[0],
                                   this->offsets[0], this->mates[0], this->seqs[0], this->quals[0],
                                   this->tag_dicts[0]),
                 seqan3::format_error);

    // no negative values except -1 are allowed fot the ref offset
    EXPECT_THROW(fout.emplace_back(&(this->header), this->ids[0], this->flags[0], this->ref_id,
                                   -3, this->mapqs[0], this->alignments[0],
                                   this->offsets[0], this->mates[0], this->seqs[0], this->quals[0],
                                   this->tag_dicts[0]),
                 seqan3::format_error);
}

REGISTER_TYPED_TEST_SUITE_P(alignment_file_read,
                            input_concept,
                            header_sucess,
                            read_in_all_data,
                            read_in_all_but_empty_data,
                            read_in_almost_nothing,
                            read_in_alignment_only_with_ref,
                            read_in_alignment_only_without_ref,
                            read_mate_but_not_ref_id_with_ref,
                            read_mate_but_not_ref_id_without_ref,
                            cigar_vector,
                            format_error_ref_id_not_in_reference_information);

REGISTER_TYPED_TEST_SUITE_P(alignment_file_write,
                            write_empty_members,
                            output_concept,
                            default_options_all_members_specified,
                            write_ref_id_with_different_types,
                            with_header,
                            cigar_vector,
                            special_cases,
                            format_errors);
