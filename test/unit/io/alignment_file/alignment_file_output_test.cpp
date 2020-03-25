// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/view/zip.hpp>
#include <range/v3/view/filter.hpp>

#include <seqan3/io/alignment_file/input.hpp>
#include <seqan3/io/alignment_file/output.hpp>
#include <seqan3/test/tmp_filename.hpp>
#include <seqan3/std/iterator>

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;

using default_fields = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::qual>;

std::vector<seqan3::dna5_vector> seqs
{
    "ACGT"_dna5,
    "AGGCTGNAGGCTGNA"_dna5,
    "GGAGTATAATATATATATATATAT"_dna5
};

std::vector<std::string> ids
{
    "read1",
    "read2",
    "read3"
};

std::string const output_comp
{
    "read1\t0\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n"
    "read2\t0\t*\t0\t0\t*\t*\t0\t0\tAGGCTGNAGGCTGNA\t*\n"
    "read3\t0\t*\t0\t0\t*\t*\t0\t0\tGGAGTATAATATATATATATATAT\t*\n"
};

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TEST(alignment_file_output_iterator, concepts)
{
    using it_t = typename seqan3::alignment_file_output<>::iterator;
    using sen_t = typename seqan3::alignment_file_output<>::sentinel;

    EXPECT_TRUE((std::output_iterator<it_t, std::tuple<std::string, std::string>>));
    EXPECT_TRUE((std::sentinel_for<sen_t, it_t>));
}

TEST(general, concepts)
{
    using t = seqan3::alignment_file_output<>;
    EXPECT_TRUE((std::ranges::output_range<t, std::tuple<std::string, std::string>>));

    using ct = seqan3::alignment_file_output<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::output_range<ct, std::tuple<std::string, std::string>>));
}

TEST(general, construct_by_filename)
{
    /* just the filename */
    {
        seqan3::test::tmp_filename filename{"alignment_file_output_constructor.sam"};
        EXPECT_NO_THROW( seqan3::alignment_file_output<>{filename.get_path()} );
    }

    /* wrong extension */
    {
        seqan3::test::tmp_filename filename{"alignment_file_output_constructor.xyz"};
        std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        EXPECT_THROW( seqan3::alignment_file_output<>{filename.get_path()},
                      seqan3::unhandled_extension_error );
    }

    /* unknown file */
    {
        seqan3::test::tmp_filename filename{"I/do/not/exist.sam"};
        EXPECT_THROW( seqan3::alignment_file_output<>{filename.get_path()}, seqan3::file_open_error );
    }

    /* filename + fields */
    using fields_seq = seqan3::fields<seqan3::field::seq>;
    {
        seqan3::test::tmp_filename filename{"alignment_file_output_constructor.sam"};
        EXPECT_NO_THROW(( seqan3::alignment_file_output<fields_seq,
                                                        seqan3::type_list<seqan3::format_sam>>{filename.get_path(),
                                                                                                fields_seq{}} ));
    }
}

TEST(general, construct_from_stream)
{
    /* stream + format_tag */
    EXPECT_NO_THROW(( seqan3::alignment_file_output<default_fields,
                                                    seqan3::type_list<seqan3::format_sam>>{std::ostringstream{},
                                                                                           seqan3::format_sam{}} ));


    /* stream + format_tag + fields */
    EXPECT_NO_THROW(( seqan3::alignment_file_output<default_fields,
                                                    seqan3::type_list<seqan3::format_sam>>{std::ostringstream{},
                                                                                           seqan3::format_sam{},
                                                                                           default_fields{}} ));
}

TEST(general, default_template_args_and_deduction_guides)
{
    using comp1 = seqan3::fields<seqan3::field::seq,
                                 seqan3::field::id,
                                 seqan3::field::offset,
                                 seqan3::field::ref_seq,
                                 seqan3::field::ref_id,
                                 seqan3::field::ref_offset,
                                 seqan3::field::alignment,
                                 seqan3::field::cigar,
                                 seqan3::field::mapq,
                                 seqan3::field::qual,
                                 seqan3::field::flag,
                                 seqan3::field::mate,
                                 seqan3::field::tags,
                                 seqan3::field::evalue,
                                 seqan3::field::bit_score,
                                 seqan3::field::header_ptr>;
    using comp2 = seqan3::type_list<seqan3::format_sam, seqan3::format_bam>;
    using comp3 = char;

    /* default template args */
    {
        using t = seqan3::alignment_file_output<>;
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided filename constructor */
    {
        seqan3::test::tmp_filename filename{"alignment_file_output_constructor.sam"};

        seqan3::alignment_file_output fout{filename.get_path()};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided filename constructor + custom fields */
    {
        seqan3::test::tmp_filename filename{"alignment_file_output_constructor.sam"};

        seqan3::alignment_file_output fout{filename.get_path(), seqan3::fields<seqan3::field::alignment>{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, seqan3::fields<seqan3::field::alignment>>)); // changed
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided stream constructor */
    {
        std::ostringstream ext{};
        seqan3::alignment_file_output fout{ext, seqan3::format_sam{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      seqan3::type_list<seqan3::format_sam>>)); // changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided stream temporary constructor */
    {
        seqan3::alignment_file_output fout{std::ostringstream{}, seqan3::format_sam{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      seqan3::type_list<seqan3::format_sam>>)); // changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }
}

// ----------------------------------------------------------------------------
// *impl
// ----------------------------------------------------------------------------

template <typename fn_t>
void row_wise_impl(fn_t fn)
{
    seqan3::alignment_file_output fout{std::ostringstream{}, seqan3::format_sam{}};

    for (size_t i = 0; i < 3; ++i)
        fn(fout, i);

    fout.get_stream().flush();
    EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), output_comp);
}

template <typename source_t>
void assign_impl(source_t && source)
{
    seqan3::alignment_file_output fout{std::ostringstream{}, seqan3::format_sam{}};

    fout = source;

    fout.get_stream().flush();
    EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), output_comp);
}

// ----------------------------------------------------------------------------
// row
// ----------------------------------------------------------------------------

TEST(row, assign_to_iterator)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>,
                       seqan3::fields<seqan3::field::seq, seqan3::field::id>> r{seqs[i], ids[i]};

        std::ranges::begin(file) = r;
    });
}

TEST(row, push_back_record)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>,
                       seqan3::fields<seqan3::field::seq, seqan3::field::id>> r{seqs[i], ids[i]};

        file.push_back(r);
    });
}

TEST(row, push_back_record_rvalue)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>,
                       seqan3::fields<seqan3::field::seq, seqan3::field::id>> r{seqs[i], ids[i]};

        file.push_back(std::move(r));
    });
}

TEST(row, push_back_record_const)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>,
                       seqan3::fields<seqan3::field::seq, seqan3::field::id>> const r{seqs[i], ids[i]};

        file.push_back(r);
    });
}

TEST(row, push_back_record_const_element)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        seqan3::record<seqan3::type_list<seqan3::dna5_vector const, std::string const>,
                       seqan3::fields<seqan3::field::seq, seqan3::field::id>> const r{seqs[i], ids[i]};

        file.push_back(r);
    });
}

TEST(row, push_back_tuple)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        std::tuple<seqan3::dna5_vector, std::string> t{seqs[i], ids[i]};

        file.push_back(t);
    });
}

TEST(row, push_back_tuple_rvalue)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        std::tuple<seqan3::dna5_vector, std::string> t{seqs[i], ids[i]};

        file.push_back(std::move(t));
    });
}

TEST(row, push_back_tuple_const)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        std::tuple<seqan3::dna5_vector, std::string> const t{seqs[i], ids[i]};

        file.push_back(t);
    });
}

TEST(row, push_back_tuple_const_element)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        std::tuple<seqan3::dna5_vector const, std::string const> t{seqs[i], ids[i]};

        file.push_back(t);
    });
}

TEST(row, emplace_back)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        file.emplace_back(seqs[i], ids[i]);
    });
}

/* Here the record contains a different field composite than the file. The record knows about the
 * association of values and fields, so it does not need to be guessed from the file.
 */
TEST(row, different_fields_in_record_and_file)
{
    std::vector<seqan3::phred42> qual;
    qual.resize(seqs[1].size());

    seqan3::record<seqan3::type_list<std::vector<seqan3::phred42>, std::string, seqan3::dna5_vector>,
                   seqan3::fields<seqan3::field::qual, seqan3::field::id, seqan3::field::seq>> rec{qual,
                                                                                                   ids[1],
                                                                                                   seqs[1]};

    seqan3::alignment_file_output fout{std::ostringstream{}, seqan3::format_sam{}, seqan3::fields<seqan3::field::seq,
                                                                                                  seqan3::field::id>{}};

    fout.emplace_back("AGGCTGNAGGCTGNA"_dna5, std::string("read1"));
    // fout.emplace_back("AGGCTGNAGGCTGNA"_dna5, "read1"); // const char * is not allowed
    fout.push_back(rec);

    fout.get_stream().flush();

    std::string const expected_out
    {
        "read1\t0\t*\t0\t0\t*\t*\t0\t0\tAGGCTGNAGGCTGNA\t*\n"
        "read2\t0\t*\t0\t0\t*\t*\t0\t0\tAGGCTGNAGGCTGNA\t!!!!!!!!!!!!!!!\n"
    };

    EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), expected_out);
}

TEST(row, print_header_in_file)
{
    std::vector<std::string> ref_ids{"ref1", "ref2"};
    std::vector<int32_t>     ref_len{234511, 243243};

    seqan3::alignment_file_output fout{std::ostringstream{},
                                       ref_ids,
                                       ref_len,
                                       seqan3::format_sam{},
                                       seqan3::fields<seqan3::field::id>{}};

    fout.emplace_back(std::string("read1"));

    fout.get_stream().flush();

    std::string const expected_out
    {
        "@HD\tVN:1.6\n"
        "@SQ\tSN:ref1\tLN:234511\n"
        "@SQ\tSN:ref2\tLN:243243\n"
        "read1\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*\n" // empty read
    };

    EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), expected_out);
}

TEST(row, print_header_in_record)
{
    std::vector<std::string> const ref_ids{"ref1", "ref2"};
    std::vector<int32_t>     const ref_len{234511, 243243};

    seqan3::alignment_file_header header{ref_ids};

    header.ref_id_info.push_back({ref_len[0], ""});
    header.ref_id_info.push_back({ref_len[1], ""});
    header.ref_dict[ref_ids[0]] = 0;
    header.ref_dict[ref_ids[1]] = 1;

    // no file header present
    {
        seqan3::alignment_file_output fout{std::ostringstream{},
                                           seqan3::format_sam{},
                                           seqan3::fields<seqan3::field::header_ptr>{}};

        fout.emplace_back(&header);

        fout.get_stream().flush();

        std::string const expected_out
        {
           "@HD\tVN:1.6\n"
           "@SQ\tSN:ref1\tLN:234511\n"
           "@SQ\tSN:ref2\tLN:243243\n"
            "*\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*\n" // empty read
        };

        EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), expected_out);
    }

    // file header present but record header pointer is favoured
    {
        seqan3::alignment_file_output fout{std::ostringstream{},
                                           std::vector<std::string>{"other_ref1", "other_ref2"},
                                           std::vector<int32_t>{12, 13},
                                           seqan3::format_sam{},
                                           seqan3::fields<seqan3::field::header_ptr>{}};

        fout.emplace_back(&header);

        fout.get_stream().flush();

        std::string const expected_out
        {
           "@HD\tVN:1.6\n"
           "@SQ\tSN:ref1\tLN:234511\n"
           "@SQ\tSN:ref2\tLN:243243\n"
            "*\t0\t*\t0\t0\t*\t*\t0\t0\t*\t*\n" // empty read
        };

        EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), expected_out);
    }
}

// ----------------------------------------------------------------------------
// rows
// ----------------------------------------------------------------------------

TEST(rows, assign_range_of_records)
{
    std::vector<seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>,
                               seqan3::fields<seqan3::field::seq, seqan3::field::id>>> range;

    for (size_t i = 0; i < 3; ++i)
        range.emplace_back(seqs[i], ids[i]);

    assign_impl(range);
}

TEST(rows, assign_range_of_records_const)
{
    std::vector<seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>,
                               seqan3::fields<seqan3::field::seq, seqan3::field::id>>> range;

    for (size_t i = 0; i < 3; ++i)
        range.emplace_back(seqs[i], ids[i]);

    assign_impl(std::as_const(range));
}

TEST(rows, assign_range_of_tuples)
{
    std::vector<std::tuple<seqan3::dna5_vector, std::string>> range;

    for (size_t i = 0; i < 3; ++i)
        range.emplace_back(seqs[i], ids[i]);

    assign_impl(range);
}

TEST(rows, assign_alignment_file_input)
{
    std::vector<std::string> ref_ids{"ref"};
    std::vector<seqan3::dna4_vector> ref_seqs{"ACTAGCTAGGAGGACTAGCATCGATC"_dna4};

    std::string comp =
R"(@HD	VN:1.6	SO:unknown	GO:none
@SQ	SN:ref	LN:26
@PG	ID:prog1	PN:cool_program
@CO	This is a comment.
read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read2	42	ref	2	62	7M1D1M1S	ref	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1D1M1I1M1I1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)";

    // with reference information
    {
        seqan3::alignment_file_input fin{std::istringstream{comp}, ref_ids, ref_seqs, seqan3::format_sam{}};
        seqan3::alignment_file_output fout{std::ostringstream{}, seqan3::format_sam{}};

        fout = fin;

        fout.get_stream().flush();

        EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), comp);
    }

    // without reference information
    {
        seqan3::alignment_file_input fin{std::istringstream{comp}, seqan3::format_sam{}};
        seqan3::alignment_file_output fout{std::ostringstream{}, seqan3::format_sam{}};

        fout = fin;

        fout.get_stream().flush();

        EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), comp);
    }
}

TEST(rows, assign_alignment_file_pipes)
{
    std::vector<std::string> const ref_ids{"ref"};
    std::vector<seqan3::dna4_vector> const ref_seqs{"ACTAGCTAGGAGGACTAGCATCGATC"_dna4};

    std::string comp =
R"(@HD	VN:1.6	SO:unknown	GO:none
@SQ	SN:ref	LN:26
@PG	ID:prog1	PN:cool_program
@CO	This is a comment.
read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read2	42	ref	2	62	7M1D1M1S	ref	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1D1M1I1M1I1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)";

    seqan3::alignment_file_input fin{std::istringstream{comp}, ref_ids, ref_seqs, seqan3::format_sam{}};
    seqan3::alignment_file_output fout{std::ostringstream{}, seqan3::format_sam{}};

    fin | fout;

    fout.get_stream().flush();

    EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), comp);
}

#if SEQAN3_HAS_ZLIB
TEST(rows, write_bam_file)
{
    seqan3::test::tmp_filename const filename{"in_out.bam"};

    std::vector<std::string> const ref_ids{"ref"};
    std::vector<seqan3::dna4_vector> const ref_seqs{"ACTAGCTAGGAGGACTAGCATCGATC"_dna4};

    std::string comp =
R"(@HD	VN:1.6	SO:unknown	GO:none
@SQ	SN:ref	LN:26
@PG	ID:prog1	PN:cool_program
@CO	This is a comment.
read1	41	ref	1	61	1S1M1D1M1I	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read2	42	ref	2	62	7M1D1M1S	ref	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1D1M1I1M1I1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)";
    {
        seqan3::alignment_file_input fin{std::istringstream{comp}, ref_ids, ref_seqs, seqan3::format_sam{}};
        seqan3::alignment_file_output fout{filename.get_path()};

        fin | fout;
    }

    seqan3::alignment_file_input fin2{filename.get_path(), ref_ids, ref_seqs};
    seqan3::alignment_file_output fout2{std::ostringstream{}, seqan3::format_sam{}};

    fin2 | fout2;

    fout2.get_stream().flush();

    EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout2.get_stream()).str(), comp);
}
#endif // SEQAN3_HAS_ZLIB

TEST(rows, convert_sam_to_blast)
{
    // TODO when blast format is implemented
}

// ----------------------------------------------------------------------------
// compression
// ----------------------------------------------------------------------------

std::string compression_by_filename_impl(seqan3::test::tmp_filename & filename)
{
    {
        // explicitly only test compression on sam format
        seqan3::alignment_file_output<typename seqan3::alignment_file_output<>::selected_field_ids,
                                      seqan3::type_list<seqan3::format_sam>,
                                      seqan3::ref_info_not_given> fout{filename.get_path()};

        for (size_t i = 0; i < 3; ++i)
        {
            seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>,
                           seqan3::fields<seqan3::field::seq, seqan3::field::id>> r{seqs[i], ids[i]};

            fout.push_back(r);
        }

    }

    std::string buffer{};

    {
        std::ifstream fi{filename.get_path(), std::ios::binary};

        buffer = std::string{std::istreambuf_iterator<char>{fi}, std::istreambuf_iterator<char>{}};
    }

    return buffer;
}

template <typename comp_stream_t>
void compression_by_stream_impl(comp_stream_t & stream)
{
    seqan3::alignment_file_output fout{stream, seqan3::format_sam{}};

    for (size_t i = 0; i < 3; ++i)
    {
        seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>,
                       seqan3::fields<seqan3::field::seq, seqan3::field::id>> r{seqs[i], ids[i]};

        fout.push_back(r);
    }
}

#ifdef SEQAN3_HAS_ZLIB
std::string expected_gz
{
    '\x1F', '\x8B', '\x08', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x2B', '\x4A', '\x4D',
    '\x4C', '\x31', '\xE4', '\x34', '\xE0', '\xD4', '\x02', '\x62', '\x10', '\x09', '\xA1', '\x1D', '\x9D',
    '\xDD', '\x43', '\x38', '\xB5', '\xB8', '\x8A', '\x80', '\x92', '\x46', '\x98', '\x92', '\xEE', '\xEE',
    '\xCE', '\x21', '\xEE', '\x7E', '\x30', '\x0A', '\xAA', '\xCE', '\x18', '\x43', '\x9D', '\xBB', '\xBB',
    '\xA3', '\x7B', '\x88', '\x63', '\x88', '\xA3', '\x63', '\x08', '\x2A', '\x04', '\x6A', '\x00', '\x00',
    '\x7E', '\x6C', '\x6C', '\x0F', '\x76', '\x00', '\x00', '\x00'
};

std::string expected_bgzf
{
    '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x06', '\x00', '\x42', '\x43',
    '\x02', '\x00', '\x50', '\x00', '\x2B', '\x4A', '\x4D', '\x4C', '\x31', '\xE4', '\x34', '\xE0', '\xD4', '\x02',
    '\x62', '\x10', '\x09', '\xA1', '\x1D', '\x9D', '\xDD', '\x43', '\x38', '\xB5', '\xB8', '\x8A', '\x80', '\x92',
    '\x46', '\x98', '\x92', '\xEE', '\xEE', '\xCE', '\x21', '\xEE', '\x7E', '\x8E', '\x50', '\x0A', '\xAA', '\xCE',
    '\x18', '\x43', '\x9D', '\xBB', '\xBB', '\xA3', '\x7B', '\x88', '\x63', '\x88', '\x23', '\x10', '\xA1', '\x40',
    '\xA0', '\x06', '\x00', '\x7E', '\x6C', '\x6C', '\x0F', '\x76', '\x00', '\x00', '\x00', '\x1F', '\x8B', '\x08',
    '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43', '\x02', '\x00', '\x1B',
    '\x00', '\x03', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00'
};

TEST(compression, by_filename_gz)
{
    seqan3::test::tmp_filename filename{"alignment_file_output_test.sam.gz"};

    std::string buffer = compression_by_filename_impl(filename);
    buffer[9] = '\x00'; // zero out OS byte.
    EXPECT_EQ(buffer, expected_bgzf);
}

TEST(compression, by_stream_gz)
{
    std::ostringstream out;

    {
        seqan3::contrib::gz_ostream compout{out};
        compression_by_stream_impl(compout);
    }
    std::string buffer = out.str();
    buffer[9] = '\x00'; // zero out OS byte.
    EXPECT_EQ(buffer, expected_gz);
}

TEST(compression, by_filename_bgzf)
{
    seqan3::test::tmp_filename filename{"alignment_file_output_test.sam.bgzf"};

    std::string buffer = compression_by_filename_impl(filename);
    buffer[9] = '\x00'; // zero out OS byte.
    EXPECT_EQ(buffer, expected_bgzf);
}

TEST(compression, by_stream_bgzf)
{
    std::ostringstream out;

    {
        seqan3::contrib::bgzf_ostream compout{out};
        compression_by_stream_impl(compout);
    }
    std::string buffer = out.str();
    buffer[9] = '\x00'; // zero out OS byte.
    EXPECT_EQ(buffer, expected_bgzf);
}
#endif

#ifdef SEQAN3_HAS_BZIP2
std::string expected_bz2
{
    '\x42', '\x5A', '\x68', '\x39', '\x31', '\x41', '\x59', '\x26', '\x53', '\x59', '\xEA', '\x2B', '\x97',
    '\x64', '\x00', '\x00', '\x39', '\xDF', '\x80', '\x00', '\x30', '\x00', '\x10', '\x78', '\x00', '\x28',
    '\x81', '\x04', '\x00', '\x26', '\x00', '\x10', '\x00', '\x20', '\x00', '\x48', '\x45', '\x4D', '\xAA',
    '\x31', '\x0C', '\x80', '\xC5', '\x19', '\x06', '\x86', '\x48', '\x31', '\xF0', '\xCC', '\x6F', '\x8C',
    '\xDC', '\x78', '\x1B', '\x38', '\x51', '\xDB', '\xAE', '\xA5', '\x5B', '\x50', '\x0E', '\xCA', '\x49',
    '\x44', '\x35', '\x4C', '\x12', '\x41', '\x20', '\x6C', '\x24', '\xC9', '\xA3', '\x47', '\xE2', '\xEE',
    '\x48', '\xA7', '\x0A', '\x12', '\x1D', '\x45', '\x72', '\xEC', '\x80'
};

TEST(compression, by_filename_bz2)
{
    seqan3::test::tmp_filename filename{"alignment_file_output_test.sam.bz2"};

    std::string buffer = compression_by_filename_impl(filename);
    EXPECT_EQ(buffer, expected_bz2);
}

TEST(compression, by_stream_bz2)
{
    std::ostringstream out;

    {
        seqan3::contrib::bz2_ostream compout{out};
        compression_by_stream_impl(compout);
    }

    EXPECT_EQ(out.str(), expected_bz2);
}
#endif
