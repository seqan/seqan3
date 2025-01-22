// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <sstream>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/debug_stream/byte.hpp>
#include <seqan3/core/debug_stream/optional.hpp>
#include <seqan3/core/debug_stream/tuple.hpp>
#include <seqan3/core/debug_stream/variant.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sam_file/input_format_concept.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/io/sam_file/output_format_concept.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/pretty_printing.hpp>
#include <seqan3/test/streambuf.hpp>

using seqan3::operator""_cigar_operation;
using seqan3::operator""_dna5;
using seqan3::operator""_phred42;
using seqan3::operator""_tag;

// global variables for reuse
seqan3::sam_file_input_options<seqan3::dna5> input_options;
seqan3::sam_file_output_options output_options;

struct sam_file_data : public ::testing::Test
{
    sam_file_data()
    {
        ref_sequences = std::vector<seqan3::dna5_vector>{ref_seq};
        ref_ids = std::vector<std::string>{ref_id};
        header = seqan3::sam_file_header{ref_ids};
        header.ref_id_info.emplace_back(ref_seq.size(), "");
        header.ref_dict[header.ref_ids()[0]] = 0; // set up header which is otherwise done on file level
    }

    // expected data for 3 reads
    // -------------------------------------------------------------------------

    std::vector<seqan3::dna5_vector> seqs{"ACGT"_dna5, "AGGCTGNAG"_dna5, "GGAGTATA"_dna5};

    std::vector<std::string> ids{"read1", "read2", "read3"};

    std::vector<std::vector<seqan3::phred42>> quals{
        {"!##$"_phred42},
        {"!##$&'()*"_phred42},
        {"!!*+,-./"_phred42},
    };

    seqan3::dna5_vector ref_seq = "ACTGATCGAGAGGATCTAGAGGAGATCGTAGGAC"_dna5;

    std::vector<seqan3::gapped<seqan3::dna5>> ref_seq_gapped1 = {'A'_dna5, 'C'_dna5, 'T'_dna5, seqan3::gap{}};
    std::vector<seqan3::gapped<seqan3::dna5>> ref_seq_gapped2 =
        {'C'_dna5, 'T'_dna5, 'G'_dna5, 'A'_dna5, 'T'_dna5, 'C'_dna5, 'G'_dna5, 'A'_dna5, 'G'_dna5};
    std::vector<seqan3::gapped<seqan3::dna5>> ref_seq_gapped3 =
        {'T'_dna5, seqan3::gap{}, 'G'_dna5, seqan3::gap{}, 'A'_dna5, seqan3::gap{}, 'T'_dna5, 'C'_dna5};

    std::vector<std::vector<seqan3::cigar>> cigars{{{1, 'S'_cigar_operation}, // 1S1M1D1M1I
                                                    {1, 'M'_cigar_operation},
                                                    {1, 'D'_cigar_operation},
                                                    {1, 'M'_cigar_operation},
                                                    {1, 'I'_cigar_operation}},

                                                   {{1, 'H'_cigar_operation}, // 1H7M1D1M1S2H
                                                    {7, 'M'_cigar_operation},
                                                    {1, 'D'_cigar_operation},
                                                    {1, 'M'_cigar_operation},
                                                    {1, 'S'_cigar_operation},
                                                    {2, 'H'_cigar_operation}},

                                                   {{1, 'S'_cigar_operation}, //1S1M1P1M1I1M1I1D1M1S
                                                    {1, 'M'_cigar_operation},
                                                    {1, 'P'_cigar_operation},
                                                    {1, 'M'_cigar_operation},
                                                    {1, 'I'_cigar_operation},
                                                    {1, 'M'_cigar_operation},
                                                    {1, 'I'_cigar_operation},
                                                    {1, 'D'_cigar_operation},
                                                    {1, 'M'_cigar_operation},
                                                    {1, 'S'_cigar_operation}}};

    std::string ref_id = "ref";

    std::vector<int32_t> ref_offsets{0, 1, 2};

    std::vector<seqan3::sam_flag> flags{seqan3::sam_flag{41u}, seqan3::sam_flag{42u}, seqan3::sam_flag{43u}};

    std::vector<uint8_t> mapqs{61u, 62u, 63u};

    std::vector<std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t>> mates{{0, 9, 300},
                                                                                           {0, 9, 300},
                                                                                           {0, 9, 300}};

    std::vector<seqan3::sam_tag_dictionary> tag_dicts = []()
    {
        std::vector<seqan3::sam_tag_dictionary> td{{}, {}, {}};
        td[0]["NM"_tag] = 7;
        td[0]["AS"_tag] = 2;
        td[1]["xy"_tag] = std::vector<uint16_t>{3, 4, 5};
        return td;
    }();

    std::vector<seqan3::sam_tag_dictionary> full_tag_dicts = []()
    {
        std::vector<seqan3::sam_tag_dictionary> td{{}, {}, {}};
        td[0]["NM"_tag] = -7;
        td[0]["AS"_tag] = 2;
        td[0]["CC"_tag] = 300;
        td[0]["cc"_tag] = -300;
        td[0]["aa"_tag] = 'c';
        td[0]["ff"_tag] = 3.1f;
        td[0]["zz"_tag] = "str";
        td[1]["bc"_tag] = std::vector<int8_t>{-3};
        td[1]["bC"_tag] = std::vector<uint8_t>{3u, 200u};
        td[1]["bs"_tag] = std::vector<int16_t>{-3, 200, -300};
        td[1]["bS"_tag] = std::vector<uint16_t>{300u, 40u, 500u};
        td[1]["bi"_tag] = std::vector<int32_t>{-3, 200, -66000};
        td[1]["bI"_tag] = std::vector<uint32_t>{294'967'296u};
        td[1]["bf"_tag] = std::vector<float>{3.5f, 0.1f, 43.8f};
        return td;
    }();

    std::vector<seqan3::dna5_vector> ref_sequences{};
    std::vector<std::string> ref_ids{};
    seqan3::sam_file_header<std::vector<std::string>> header{};
};

template <typename format_t>
struct sam_file_read : public sam_file_data
{};

TYPED_TEST_SUITE_P(sam_file_read);

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TYPED_TEST_P(sam_file_read, input_concept)
{
    EXPECT_TRUE((seqan3::sam_file_input_format<TypeParam>));
}

// ----------------------------------------------------------------------------
// sam_file_read
// ----------------------------------------------------------------------------

TYPED_TEST_P(sam_file_read, header_sucess)
{
    typename TestFixture::stream_type istream{this->big_header_input};

    seqan3::sam_file_input fin{istream, TypeParam{}};
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

    EXPECT_EQ(header.ref_id_info[(header.ref_dict[id1])], (std::tuple<uint32_t, std::string>{249'250'621u, ""}));
    EXPECT_EQ(header.ref_id_info[(header.ref_dict[id2])],
              (std::tuple<uint32_t, std::string>{243'199'373u, "AS:hs37d5"}));

    EXPECT_EQ(header.read_groups[0],
              (std::pair<std::string, std::string>{"U0a_A2_L1", "PL:illumina\tPU:1\tLB:1\tSM:NA12878"}));
    EXPECT_EQ(header.read_groups[1],
              (std::pair<std::string, std::string>{"U0a_A2_L2", "PL:illumina\tSM:NA12878\tPU:1\tLB:1"}));

    EXPECT_EQ(header.comments[0], "Tralalalalalala this is a comment");
}

TYPED_TEST_P(sam_file_read, read_in_all_data)
{
    typename TestFixture::stream_type istream{this->verbose_reads_input};
    seqan3::sam_file_input fin{istream, this->ref_ids, this->ref_sequences, TypeParam{}};

    this->full_tag_dicts[1]["bH"_tag] = std::vector<std::byte>{std::byte{0x1A}, std::byte{0xE3}, std::byte{0x01}};

    size_t i{0};
    for (auto & rec : fin)
    {
        EXPECT_EQ(rec.sequence(), this->seqs[i]);
        EXPECT_EQ(rec.id(), this->ids[i]);
        EXPECT_EQ(rec.base_qualities(), this->quals[i]);
        EXPECT_RANGE_EQ(rec.cigar_sequence(), this->cigars[i]);
        EXPECT_EQ(rec.reference_id(), 0);
        EXPECT_EQ(*rec.reference_position(), this->ref_offsets[i]);
        EXPECT_EQ(rec.flag(), this->flags[i]);
        EXPECT_EQ(rec.mapping_quality(), this->mapqs[i]);
        EXPECT_EQ(rec.mate_reference_id(), std::get<0>(this->mates[i]));
        EXPECT_EQ(rec.mate_position(), std::get<1>(this->mates[i]));
        EXPECT_EQ(rec.template_length(), std::get<2>(this->mates[i]));
        EXPECT_EQ(rec.tags(), this->full_tag_dicts[i]);
        ++i;
    }
}

TYPED_TEST_P(sam_file_read, read_in_all_data_with_small_stream_buffer)
{
    typename TestFixture::stream_type istream{this->verbose_reads_input};
    std::istream & in{istream};
    std::streambuf * orig = in.rdbuf();
    seqan3::test::streambuf_with_custom_buffer_size<20> buf(orig);
    in.rdbuf(&buf);

    seqan3::sam_file_input fin{in, this->ref_ids, this->ref_sequences, TypeParam{}};

    this->full_tag_dicts[1]["bH"_tag] = std::vector<std::byte>{std::byte{0x1A}, std::byte{0xE3}, std::byte{0x01}};

    size_t i{0};
    for (auto & rec : fin)
    {
        EXPECT_EQ(rec.sequence(), this->seqs[i]);
        EXPECT_EQ(rec.id(), this->ids[i]);
        EXPECT_EQ(rec.base_qualities(), this->quals[i]);
        EXPECT_RANGE_EQ(rec.cigar_sequence(), this->cigars[i]);
        EXPECT_EQ(rec.reference_id(), 0);
        EXPECT_EQ(*rec.reference_position(), this->ref_offsets[i]);
        EXPECT_EQ(rec.flag(), this->flags[i]);
        EXPECT_EQ(rec.mapping_quality(), this->mapqs[i]);
        EXPECT_EQ(rec.mate_reference_id(), std::get<0>(this->mates[i]));
        EXPECT_EQ(rec.mate_position(), std::get<1>(this->mates[i]));
        EXPECT_EQ(rec.template_length(), std::get<2>(this->mates[i]));
        EXPECT_EQ(rec.tags(), this->full_tag_dicts[i]);
        ++i;
    }
}

TYPED_TEST_P(sam_file_read, read_in_all_but_empty_data)
{
    typename TestFixture::stream_type istream{this->empty_input};
    seqan3::sam_file_input fin{istream, this->ref_ids, this->ref_sequences, TypeParam{}};

    EXPECT_TRUE((*fin.begin()).sequence().empty());
    EXPECT_TRUE((*fin.begin()).id().empty());
    EXPECT_TRUE((*fin.begin()).base_qualities().empty());
    EXPECT_TRUE((*fin.begin()).cigar_sequence().empty());
    EXPECT_TRUE(!(*fin.begin()).reference_id().has_value());
    EXPECT_TRUE(!(*fin.begin()).reference_position().has_value());
    EXPECT_EQ((*fin.begin()).flag(), seqan3::sam_flag{0u});
    EXPECT_EQ((*fin.begin()).mapping_quality(), 0u);
    EXPECT_TRUE(!(*fin.begin()).mate_reference_id().has_value());
    EXPECT_TRUE(!(*fin.begin()).mate_position().has_value());
    EXPECT_EQ((*fin.begin()).template_length(), int32_t{});
    EXPECT_TRUE((*fin.begin()).tags().empty());
}

TYPED_TEST_P(sam_file_read, read_in_almost_nothing)
{
    typename TestFixture::stream_type istream{this->simple_three_reads_input};
    seqan3::sam_file_input fin{istream, TypeParam{}, seqan3::fields<seqan3::field::mapq>{}};

    size_t i{0};
    for (auto & [mapq] : fin)
        EXPECT_EQ(mapq, this->mapqs[i++]);
}

TYPED_TEST_P(sam_file_read, read_mate_but_not_ref_id_with_ref)
{
#if __cplusplus > 202002L && defined(__clang__)
    GTEST_SKIP() << "Weird error with clang and CPP23 tuples";
#else
    { /*with reference information*/
        typename TestFixture::stream_type istream{this->simple_three_reads_input};
        seqan3::sam_file_input fin{istream,
                                   this->ref_ids,
                                   this->ref_sequences,
                                   TypeParam{},
                                   seqan3::fields<seqan3::field::mate>{}};

        size_t i{0};
        for (auto & [mate] : fin)
            EXPECT_EQ(mate, this->mates[i++]);
    }
#endif
}

TYPED_TEST_P(sam_file_read, read_mate_but_not_ref_id_without_ref)
{
#if __cplusplus > 202002L && defined(__clang__)
    GTEST_SKIP() << "Weird error with clang and CPP23 tuples";
#else
    [[maybe_unused]] std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> mate;

    { /*no reference information*/
        typename TestFixture::stream_type istream{this->simple_three_reads_input};
        seqan3::sam_file_input fin{istream, TypeParam{}, seqan3::fields<seqan3::field::mate>{}};

        size_t i{0};
        for (auto & [mate] : fin)
            EXPECT_EQ(mate, this->mates[i++]);
    }
#endif
}

TYPED_TEST_P(sam_file_read, cigar_vector)
{
    typename TestFixture::stream_type istream{this->simple_three_reads_input};
    seqan3::sam_file_input fin{istream, TypeParam{}, seqan3::fields<seqan3::field::cigar>{}};

    size_t i{0};
    for (auto & [cigar_v] : fin)
        EXPECT_EQ(cigar_v, this->cigars[i++]);
}

TYPED_TEST_P(sam_file_read, format_error_ref_id_not_in_reference_information)
{
    { // with reference information given
        typename TestFixture::stream_type istream{this->unknown_ref};
        seqan3::sam_file_input fin{istream, this->ref_ids, this->ref_sequences, TypeParam{}};
        EXPECT_THROW((fin.begin()), seqan3::format_error);
    }

    { // with reference information in the header
        typename TestFixture::stream_type istream{this->unknown_ref_header};
        seqan3::sam_file_input fin{istream, TypeParam{}};
        EXPECT_THROW((fin.begin()), seqan3::format_error);
    }
}

TYPED_TEST_P(sam_file_read, format_error_uneven_hexadecimal_tag)
{
    typename TestFixture::stream_type istream{this->wrong_hexadecimal_tag};
    seqan3::sam_file_input fin{istream, TypeParam{}};
    EXPECT_THROW((fin.begin()), seqan3::format_error);
}

// https://github.com/seqan/seqan3/pull/2423
TYPED_TEST_P(sam_file_read, issue2423)
{
    typename TestFixture::stream_type istream{this->many_refs};
    seqan3::sam_file_input fin{istream, TypeParam{}};
    ASSERT_NO_THROW(fin.begin());

    EXPECT_EQ(fin.header().ref_id_info.size(), 64u);
    EXPECT_EQ(fin.header().ref_dict.size(), 64u);
}

TYPED_TEST_P(sam_file_read, unknown_header_tag)
{
    typename TestFixture::stream_type istream{this->unknown_tag_header};
    seqan3::sam_file_input fin{istream, TypeParam{}};
    ASSERT_NO_THROW(fin.begin());

    EXPECT_EQ(fin.header().user_tags, "pb:5.0.0\totter");                        // HD
    EXPECT_EQ(std::get<1>(fin.header().ref_id_info.front()), "pb:5.0.0\totter"); // SQ
    EXPECT_EQ(std::get<1>(fin.header().read_groups.front()), "pb:5.0.0\totter"); // RG
    EXPECT_EQ(fin.header().program_infos.front().user_tags, "pb:5.0.0\totter");  // PG
}

// ----------------------------------------------------------------------------
// sam_file_write
// ----------------------------------------------------------------------------

using sam_fields = seqan3::fields<seqan3::field::header_ptr,
                                  seqan3::field::id,
                                  seqan3::field::flag,
                                  seqan3::field::ref_id,
                                  seqan3::field::ref_offset,
                                  seqan3::field::mapq,
                                  seqan3::field::cigar,
                                  seqan3::field::mate,
                                  seqan3::field::seq,
                                  seqan3::field::qual,
                                  seqan3::field::tags>;

template <typename format_type>
struct sam_file_write : public sam_file_read<format_type>
{
    std::ostringstream ostream;
};

TYPED_TEST_SUITE_P(sam_file_write);

TYPED_TEST_P(sam_file_write, output_concept)
{
    EXPECT_TRUE((seqan3::sam_file_output_format<TypeParam>));
}

TYPED_TEST_P(sam_file_write, no_records)
{
    {
        auto ref_lengths = this->ref_sequences
                         | std::views::transform(
                               [](auto const & v)
                               {
                                   return v.size();
                               });
        seqan3::sam_file_output fout{this->ostream, this->ref_ids, ref_lengths, TypeParam{}, sam_fields{}};
    }

    this->ostream.flush();
    EXPECT_EQ(this->ostream.str(), this->minimal_header);
}

TYPED_TEST_P(sam_file_write, write_empty_members)
{
    {
        seqan3::sam_file_output fout{this->ostream, TypeParam{}, sam_fields{}};

        using default_mate_t = std::tuple<std::string_view, std::optional<int32_t>, int32_t>;

        fout.emplace_back(&(this->header),
                          std::string_view{},
                          seqan3::sam_flag::none,
                          std::string_view{},
                          -1,
                          0,
                          std::vector<seqan3::cigar>{},
                          default_mate_t{},
                          std::string_view{},
                          std::string_view{},
                          seqan3::sam_tag_dictionary{});
    }

    this->ostream.flush();
    EXPECT_EQ(this->ostream.str(), this->empty_input);
}

TYPED_TEST_P(sam_file_write, default_options_all_members_specified)
{
    {
        seqan3::sam_file_output fout{this->ostream, TypeParam{}, sam_fields{}};

        for (size_t i = 0; i < 3; ++i)
        {
            ASSERT_NO_THROW(fout.emplace_back(&(this->header),
                                              this->ids[i],
                                              this->flags[i],
                                              0 /*ref_id*/,
                                              this->ref_offsets[i],
                                              this->mapqs[i],
                                              this->cigars[i],
                                              this->mates[i],
                                              this->seqs[i],
                                              this->quals[i],
                                              this->tag_dicts[i]));
        }
    }
    this->ostream.flush();

    EXPECT_EQ(this->ostream.str(), this->simple_three_reads_input);
}

TYPED_TEST_P(sam_file_write, write_ref_id_with_different_types)
{
    {
        // header ref_id_type is std::string
        seqan3::sam_file_output fout{this->ostream, TypeParam{}, sam_fields{}};

        // std::string
        ASSERT_NO_THROW(fout.emplace_back(&(this->header),
                                          this->ids[0],
                                          this->flags[0],
                                          /*----------------------->*/ this->ref_id,
                                          this->ref_offsets[0],
                                          this->mapqs[0],
                                          this->cigars[0],
                                          this->mates[0],
                                          this->seqs[0],
                                          this->quals[0],
                                          this->tag_dicts[0]));

        // std::string_view
        ASSERT_NO_THROW(fout.emplace_back(&(this->header),
                                          this->ids[1],
                                          this->flags[1],
                                          /*----------------------->*/ std::string_view{this->ref_id},
                                          this->ref_offsets[1],
                                          this->mapqs[1],
                                          this->cigars[1],
                                          this->mates[1],
                                          this->seqs[1],
                                          this->quals[1],
                                          this->tag_dicts[1]));

        // view on string
        ASSERT_NO_THROW(fout.emplace_back(&(this->header),
                                          this->ids[2],
                                          this->flags[2],
                                          /*----------------------->*/ this->ref_id | std::views::take(20),
                                          this->ref_offsets[2],
                                          this->mapqs[2],
                                          this->cigars[2],
                                          this->mates[2],
                                          this->seqs[2],
                                          this->quals[2],
                                          this->tag_dicts[2]));
    }

    this->ostream.flush();

    EXPECT_EQ(this->ostream.str(), this->simple_three_reads_input);
}

TYPED_TEST_P(sam_file_write, with_header)
{
    seqan3::sam_file_header header{std::vector<std::string>{this->ref_id}};
    header.sorting = "unknown";
    header.grouping = "none";
    header.ref_id_info.push_back({this->ref_seq.size(), "AN:other_name\tpb:5.0.0\totter"});
    header.ref_dict[this->ref_id] = 0;
    header.program_infos.push_back({"prog1", "cool_program", "./prog1", "a", "b", "c", "pb:5.0.0\totter"});
    header.read_groups.emplace_back("group1", "DS:more info\tpb:5.0.0\totter");
    header.comments.push_back("This is a comment.");
    header.user_tags = "pb:5.0.0\totter";

    {
        seqan3::sam_file_output fout{this->ostream, TypeParam{}, sam_fields{}};

        for (size_t i = 0; i < 3; ++i)
        {
            ASSERT_NO_THROW(fout.emplace_back(&header,
                                              this->ids[i],
                                              this->flags[i],
                                              0 /*ref_id*/,
                                              this->ref_offsets[i],
                                              this->mapqs[i],
                                              this->cigars[i],
                                              this->mates[i],
                                              this->seqs[i],
                                              this->quals[i],
                                              this->full_tag_dicts[i]));
        }
    }

    this->ostream.flush();

    EXPECT_EQ(this->ostream.str(), this->verbose_output);
}

TYPED_TEST_P(sam_file_write, cigar_vector)
{
    {
        seqan3::sam_file_output fout{this->ostream, TypeParam{}};
        for (size_t i = 0ul; i < 3ul; ++i)
        {
            ASSERT_NO_THROW(fout.emplace_back(this->seqs[i],
                                              this->ids[i],
                                              0 /*ref_id*/,
                                              this->ref_offsets[i],
                                              this->cigars[i],
                                              this->mapqs[i],
                                              this->quals[i],
                                              this->flags[i],
                                              this->mates[i],
                                              this->tag_dicts[i],
                                              &(this->header)));
        }
    }

    this->ostream.flush();
    // compare to original input because hard clipping is preserved when writing the cigar vector directly
    EXPECT_EQ(this->ostream.str(), this->simple_three_reads_input);

    this->ostream = std::ostringstream{}; // clear

    // 2. Write only the cigar, not the alignment
    {
        seqan3::sam_file_output fout{this->ostream,
                                     TypeParam{},
                                     seqan3::fields<seqan3::field::header_ptr,
                                                    seqan3::field::id,
                                                    seqan3::field::flag,
                                                    seqan3::field::ref_id,
                                                    seqan3::field::ref_offset,
                                                    seqan3::field::mapq,
                                                    seqan3::field::cigar,
                                                    seqan3::field::mate,
                                                    seqan3::field::seq,
                                                    seqan3::field::qual,
                                                    seqan3::field::tags>{}};

        for (size_t i = 0ul; i < 3ul; ++i)
        {
            ASSERT_NO_THROW(fout.emplace_back(&(this->header),
                                              this->ids[i],
                                              this->flags[i],
                                              0 /*ref_id*/,
                                              this->ref_offsets[i],
                                              this->mapqs[i],
                                              this->cigars[i],
                                              this->mates[i],
                                              this->seqs[i],
                                              this->quals[i],
                                              this->tag_dicts[i]));
        }
    }

    this->ostream.flush();
    // compare to original input because hard clipping is preserved when writing the cigar vector directly
    EXPECT_EQ(this->ostream.str(), this->simple_three_reads_input);
}

TYPED_TEST_P(sam_file_write, special_cases)
{
    std::optional<int32_t> rid;

    // write an empty std::optional for ref id and mate
    {
        std::tuple<std::optional<int32_t>, std::optional<int32_t>, int32_t> mate{rid, rid, 0};

        seqan3::sam_file_output fout{this->ostream, TypeParam{}, sam_fields{}};

        // std::string
        ASSERT_NO_THROW(fout.emplace_back(&(this->header),
                                          this->ids[0],
                                          this->flags[0],
                                          rid,
                                          this->ref_offsets[0],
                                          this->mapqs[0],
                                          this->cigars[0],
                                          mate,
                                          this->seqs[0],
                                          this->quals[0],
                                          seqan3::sam_tag_dictionary{}));
    }

    this->ostream.flush();
    EXPECT_EQ(this->ostream.str(), this->special_output);

    this->ostream = std::ostringstream{}; // clear

    {
        // write the ref id and mate ref as string
        std::tuple<std::string, std::optional<int32_t>, int32_t> mate_str{"", rid, 0};

        seqan3::sam_file_output fout{this->ostream, TypeParam{}, sam_fields{}};

        // std::string
        ASSERT_NO_THROW(fout.emplace_back(&(this->header),
                                          this->ids[0],
                                          this->flags[0],
                                          std::string(""),
                                          this->ref_offsets[0],
                                          this->mapqs[0],
                                          this->cigars[0],
                                          mate_str,
                                          this->seqs[0],
                                          this->quals[0],
                                          seqan3::sam_tag_dictionary{}));
    }

    this->ostream.flush();
    EXPECT_EQ(this->ostream.str(), this->special_output);
}

TYPED_TEST_P(sam_file_write, format_errors)
{
    auto ref_lengths = this->ref_sequences
                     | std::views::transform(
                           [](auto const & v)
                           {
                               return v.size();
                           });
    seqan3::sam_file_output fout{this->ostream, this->ref_ids, ref_lengths, TypeParam{}, sam_fields{}};

    // ensure that only a ref_id that is listed in the header is allowed
    EXPECT_THROW(fout.emplace_back(&(this->header),
                                   this->ids[0],
                                   this->flags[0],
                                   std::string("ref_id_that_does_not_exist"),
                                   this->ref_offsets[0],
                                   this->mapqs[0],
                                   this->cigars[0],
                                   this->mates[0],
                                   this->seqs[0],
                                   this->quals[0],
                                   this->tag_dicts[0]),
                 seqan3::format_error);

    // no negative values except -1 are allowed for the ref offset
    EXPECT_THROW(fout.emplace_back(&(this->header),
                                   this->ids[0],
                                   this->flags[0],
                                   this->ref_id,
                                   -3,
                                   this->mapqs[0],
                                   this->cigars[0],
                                   this->mates[0],
                                   this->seqs[0],
                                   this->quals[0],
                                   this->tag_dicts[0]),
                 seqan3::format_error);
}

TYPED_TEST_P(sam_file_write, issue3299)
{
    using sam_file_output_t = seqan3::sam_file_output<typename seqan3::sam_file_output<>::selected_field_ids,
                                                      seqan3::type_list<TypeParam>,
                                                      std::vector<std::string>>;
    std::vector<std::string> seq_names{"hello", "world"};
    std::vector<size_t> seq_lengths{1000, 2000};

    // Issue: Moved-from sam_file_output would try to write header on destruction
    {
        sam_file_output_t fout1{std::ostringstream{}, seq_names, seq_lengths, TypeParam{}};
        sam_file_output_t fout2{std::move(fout1)};
    }

    // Issue: Header does not own ref_ids: ref_ids outlives sam_file_output
    {
        std::vector<sam_file_output_t> alignment_streams;
        auto seq_names_copy = seq_names;
        alignment_streams.emplace_back(std::ostringstream{}, seq_names_copy, seq_lengths, TypeParam{});
        // Destructor calls:
        // 1) seq_names_copy
        // 2) alignment_streams, starting with the one element it holds
    }

    // Issue: Header does not own ref_ids: ref_ids may change
    size_t const iterations{this->issue3299_output.size()};
    // Order of destruction of vector elements differs between GCC and Clang
    std::vector<std::ostringstream> outputs(iterations);
    {
        std::vector<sam_file_output_t> alignment_streams;
        for (size_t i = 0; i < iterations; ++i)
        {
            alignment_streams.emplace_back(outputs[i], seq_names, seq_lengths, TypeParam{});

            std::ranges::for_each(seq_names,
                                  [](std::string & str)
                                  {
                                      str += "foo";
                                  });
            std::ranges::for_each(seq_lengths,
                                  [](size_t & len)
                                  {
                                      ++len;
                                  });
        }
    }
    for (size_t i = 0; i < iterations; ++i)
    {
        EXPECT_EQ(outputs[i].str(), this->issue3299_output[i]) << "Iteration: " << i;
    }
}

REGISTER_TYPED_TEST_SUITE_P(sam_file_read,
                            input_concept,
                            header_sucess,
                            read_in_all_data,
                            read_in_all_data_with_small_stream_buffer,
                            read_in_all_but_empty_data,
                            read_in_almost_nothing,
                            read_mate_but_not_ref_id_with_ref,
                            read_mate_but_not_ref_id_without_ref,
                            cigar_vector,
                            format_error_ref_id_not_in_reference_information,
                            format_error_uneven_hexadecimal_tag,
                            issue2423,
                            unknown_header_tag);

REGISTER_TYPED_TEST_SUITE_P(sam_file_write,
                            no_records,
                            write_empty_members,
                            output_concept,
                            default_options_all_members_specified,
                            write_ref_id_with_different_types,
                            with_header,
                            cigar_vector,
                            special_cases,
                            format_errors,
                            issue3299);
