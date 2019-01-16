// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/view/zip.hpp>
#include <range/v3/view/filter.hpp>

#include <seqan3/io/alignment_file/output.hpp>
// #include <seqan3/io/alignment_file/input.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/test/tmp_filename.hpp>
#include <seqan3/std/iterator>

using namespace seqan3;

std::vector<dna5_vector> seqs
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
    using it_t = typename alignment_file_output<>::iterator;
    using sen_t = typename alignment_file_output<>::sentinel;

    EXPECT_TRUE((std::OutputIterator<it_t, std::tuple<std::string, std::string>>));
    EXPECT_TRUE((std::Sentinel<sen_t, it_t>));
}

TEST(general, concepts)
{
    using t = alignment_file_output<>;
    EXPECT_TRUE((std::ranges::OutputRange<t, std::tuple<std::string, std::string>>));

    using ct = alignment_file_output<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::OutputRange<ct, std::tuple<std::string, std::string>>));
}

TEST(general, construct_by_filename)
{
    /* just the filename */
    {
        test::tmp_filename filename{"alignment_file_output_constructor.sam"};
        EXPECT_NO_THROW( alignment_file_output<>{filename.get_path()} );
    }

    /* wrong extension */
    {
        test::tmp_filename filename{"alignment_file_output_constructor.xyz"};
        std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        EXPECT_THROW( alignment_file_output<>{filename.get_path()} ,
                      unhandled_extension_error );
    }

    /* filename + fields */
    {
        test::tmp_filename filename{"alignment_file_output_constructor.sam"};
        EXPECT_NO_THROW(( alignment_file_output<fields<field::SEQ>,
                                                type_list<alignment_file_format_sam>,
                                                std::ofstream>{filename.get_path(), fields<field::SEQ>{}} ));
    }
}

TEST(general, construct_from_stream)
{
    /* stream + format_tag */
    EXPECT_NO_THROW(( alignment_file_output<fields<field::SEQ, field::ID, field::QUAL>,
                                            type_list<alignment_file_format_sam>,
                                            std::ostringstream>{std::ostringstream{},
                                                                alignment_file_format_sam{}} ));


    /* stream + format_tag + fields */
    EXPECT_NO_THROW(( alignment_file_output<fields<field::SEQ, field::ID, field::QUAL>,
                                            type_list<alignment_file_format_sam>,
                                            std::ostringstream>{std::ostringstream{},
                                                                alignment_file_format_sam{},
                                                                fields<field::SEQ, field::ID, field::QUAL>{}} ));
}

TEST(general, default_template_args_and_deduction_guides)
{
    using comp1 = fields<field::SEQ,
                         field::ID,
                         field::OFFSET,
                         field::REF_SEQ,
                         field::REF_ID,
                         field::REF_OFFSET,
                         field::ALIGNMENT,
                         field::MAPQ,
                         field::QUAL,
                         field::FLAG,
                         field::MATE,
                         field::TAGS,
                         field::EVALUE,
                         field::BIT_SCORE,
                         field::HEADER_PTR>;
    using comp2 = type_list<alignment_file_format_sam>;
    using comp3 = std::ofstream;

    /* default template args */
    {
        using t = alignment_file_output<>;
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
    }

    /* guided filename constructor */
    {
        test::tmp_filename filename{"alignment_file_output_constructor.sam"};

        alignment_file_output fout{filename.get_path()};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
    }

    /* guided filename constructor + custom fields */
    {
        test::tmp_filename filename{"alignment_file_output_constructor.sam"};

        alignment_file_output fout{filename.get_path(), fields<field::ALIGNMENT>{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::ALIGNMENT>>));             // changed
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
    }

    /* guided stream constructor */
    {
        alignment_file_output fout{std::ostringstream{}, alignment_file_format_sam{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<alignment_file_format_sam>>)); // changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        std::ostringstream>));                   // changed
    }

    /* guided stream constructor + custom fields */
    {
        alignment_file_output fout{std::ostringstream{}, alignment_file_format_sam{}, fields<field::REF_ID>{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::REF_ID>>));                   // changed
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<alignment_file_format_sam>>));// changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        std::ostringstream>));                   // changed
    }
}

// ----------------------------------------------------------------------------
// *impl
// ----------------------------------------------------------------------------

template <typename fn_t>
void row_wise_impl(fn_t fn)
{
    alignment_file_output fout{std::ostringstream{}, alignment_file_format_sam{}};

    for (size_t i = 0; i < 3; ++i)
        fn(fout, i);

    fout.get_stream().flush();
    EXPECT_EQ(fout.get_stream().str(), output_comp);
}

template <typename source_t>
void assign_impl(source_t && source)
{
    alignment_file_output fout{std::ostringstream{}, alignment_file_format_sam{}};

    fout = source;

    fout.get_stream().flush();
    EXPECT_EQ(fout.get_stream().str(), output_comp);
}

// ----------------------------------------------------------------------------
// row
// ----------------------------------------------------------------------------

TEST(row, assign_to_iterator)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        record<type_list<dna5_vector, std::string>, fields<field::SEQ, field::ID>> r{seqs[i], ids[i]};

        begin(file) = r;
    });
}

TEST(row, push_back_record)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        record<type_list<dna5_vector, std::string>, fields<field::SEQ, field::ID>> r{seqs[i], ids[i]};

        file.push_back(r);
    });
}

TEST(row, push_back_record_rvalue)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        record<type_list<dna5_vector, std::string>, fields<field::SEQ, field::ID>> r{seqs[i], ids[i]};

        file.push_back(std::move(r));
    });
}

TEST(row, push_back_record_const)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        record<type_list<dna5_vector, std::string>, fields<field::SEQ, field::ID>> const r{seqs[i], ids[i]};

        file.push_back(r);
    });
}

TEST(row, push_back_record_const_element)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        record<type_list<dna5_vector const, std::string const>, fields<field::SEQ, field::ID>> const r{seqs[i], ids[i]};

        file.push_back(r);
    });
}

TEST(row, push_back_tuple)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        std::tuple<dna5_vector, std::string> t{seqs[i], ids[i]};

        file.push_back(t);
    });
}

TEST(row, push_back_tuple_rvalue)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        std::tuple<dna5_vector, std::string> t{seqs[i], ids[i]};

        file.push_back(std::move(t));
    });
}

TEST(row, push_back_tuple_const)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        std::tuple<dna5_vector, std::string> const t{seqs[i], ids[i]};

        file.push_back(t);
    });
}

TEST(row, push_back_tuple_const_element)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        std::tuple<dna5_vector const, std::string const> t{seqs[i], ids[i]};

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

/* Here the record contains a different field composition than the file. The record knows about the
 * association of values and fields, so it does not need to be guessed from the file.
 */
TEST(row, different_fields_in_record_and_file)
{
    std::vector<phred42> qual;
    qual.resize(seqs[1].size());

    record<type_list<std::vector<phred42>, std::string, dna5_vector>,
           fields<field::QUAL, field::ID, field::SEQ>> rec{qual, ids[1], seqs[1]};

    alignment_file_output fout{std::ostringstream{}, alignment_file_format_sam{}, fields<field::SEQ, field::ID>{}};

    fout.emplace_back("AGGCTGNAGGCTGNA"_dna5, std::string("read1"));
    // fout.emplace_back("AGGCTGNAGGCTGNA"_dna5, "read1"); // const char * is not allowed
    fout.push_back(rec);

    fout.get_stream().flush();

    std::string const expected_out
    {
        "read1\t0\t*\t0\t0\t*\t*\t0\t0\tAGGCTGNAGGCTGNA\t*\n"
        "read2\t0\t*\t0\t0\t*\t*\t0\t0\tAGGCTGNAGGCTGNA\t!!!!!!!!!!!!!!!\n"
    };

    EXPECT_EQ(fout.get_stream().str(), expected_out);
}

// ----------------------------------------------------------------------------
// rows
// ----------------------------------------------------------------------------

TEST(rows, assign_range_of_records)
{
    std::vector<record<type_list<dna5_vector, std::string>, fields<field::SEQ, field::ID>>> range;

    for (size_t i = 0; i < 3; ++i)
        range.emplace_back(seqs[i], ids[i]);

    assign_impl(range);
}

TEST(rows, assign_range_of_records_const)
{
    std::vector<record<type_list<dna5_vector, std::string>, fields<field::SEQ, field::ID>>> range;

    for (size_t i = 0; i < 3; ++i)
        range.emplace_back(seqs[i], ids[i]);

    assign_impl(std::as_const(range));
}

TEST(rows, assign_range_of_tuples)
{
    std::vector<std::tuple<dna5_vector, std::string>> range;

    for (size_t i = 0; i < 3; ++i)
        range.emplace_back(seqs[i], ids[i]);

    assign_impl(range);
}

TEST(rows, assign_alignment_file_in)
{
    // TODO when input is implemented
}

TEST(rows, assign_alignment_file_pipes)
{
    // TODO when input is implemented
}

TEST(rows, convert_sam_to_blast)
{
    // TODO when blast format is implemented
}
