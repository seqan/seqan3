// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/view/zip.hpp>
#include <range/v3/view/filter.hpp>

#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/test/tmp_filename.hpp>
#include <seqan3/std/iterator>

using namespace seqan3;

TEST(sequence_file_output_iterator, concepts)
{
    using it_t = typename sequence_file_output<>::iterator;
    using sen_t = typename sequence_file_output<>::sentinel;

    EXPECT_TRUE((std::OutputIterator<it_t, std::tuple<std::string, std::string>>));
    EXPECT_TRUE((std::Sentinel<sen_t, it_t>));
}

std::vector<dna5_vector> seqs
{
    "ACGT"_dna5,
    "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN"_dna5,
    "GGAGTATAATATATATATATATAT"_dna5
};

std::vector<std::string> ids
{
    "TEST 1",
    "Test2",
    "Test3"
};

std::string const output_comp
{
    "> TEST 1\n"
    "ACGT\n"
    "> Test2\n"
    "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN\n"
    "> Test3\n"
    "GGAGTATAATATATATATATATAT\n"
};

std::vector<std::vector<phred42>> quals
{
    "!!!!"_phred42,
    "!#@$!#@$!#@#!$@#!$@#!$!#@$!#@#!$@#!$!#$@!!$$$$$$$$$$$$!!!!!!!!!!!!!!!!!!!!$$$$$$$$$$!!!!!$!"_phred42,
    "!@#!@#!#!######@$!#@!!!@"_phred42
};

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TEST(general, concepts)
{
    using t = sequence_file_output<>;
    EXPECT_TRUE((std::ranges::OutputRange<t, std::tuple<std::string, std::string>>));

    using ct = sequence_file_output<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::OutputRange<ct, std::tuple<std::string, std::string>>));
}

TEST(general, construct_by_filename)
{
    /* just the filename */
    {
        test::tmp_filename filename{"sequence_file_output_constructor.fasta"};
        EXPECT_NO_THROW( sequence_file_output<>{filename.get_path()} );
    }

    /* wrong extension */
    {
        test::tmp_filename filename{"sequence_file_output_constructor.xyz"};
        std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        EXPECT_THROW( sequence_file_output<>{filename.get_path()} ,
                      unhandled_extension_error );
    }

    /* unknown file */
    {
        test::tmp_filename filename{"I/do/not/exist.fasta"};
        EXPECT_THROW( sequence_file_output<>{filename.get_path()}, file_open_error );
    }

    /* filename + fields */
    {
        test::tmp_filename filename{"sequence_file_output_constructor.fasta"};
        EXPECT_NO_THROW(( sequence_file_output<fields<field::SEQ>,
                                            type_list<sequence_file_format_fasta>>{filename.get_path(), fields<field::SEQ>{}} ));
    }
}

TEST(general, construct_from_stream)
{
    /* stream + format_tag */
    EXPECT_NO_THROW(( sequence_file_output<fields<field::SEQ, field::ID, field::QUAL>,
                                           type_list<sequence_file_format_fasta>>{std::ostringstream{},
                                                            sequence_file_format_fasta{}} ));


    /* stream + format_tag + fields */
    EXPECT_NO_THROW(( sequence_file_output<fields<field::SEQ, field::ID, field::QUAL>,
                                        type_list<sequence_file_format_fasta>>{std::ostringstream{},
                                                            sequence_file_format_fasta{},
                                                            fields<field::SEQ, field::ID, field::QUAL>{}} ));
}

TEST(general, default_template_args_and_deduction_guides)
{
    using comp1 = fields<field::SEQ, field::ID, field::QUAL>;
    using comp2 = type_list<sequence_file_format_embl, sequence_file_format_fasta, sequence_file_format_fastq,
                            sequence_file_format_sam>;
    using comp3 = char;

    /* default template args */
    {
        using t = sequence_file_output<>;
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided filename constructor */
    {
        test::tmp_filename filename{"sequence_file_output_constructor.fasta"};

        sequence_file_output fout{filename.get_path()};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided filename constructor + custom fields */
    {
        test::tmp_filename filename{"sequence_file_output_constructor.fasta"};

        sequence_file_output fout{filename.get_path(), fields<field::SEQ>{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));                   // changed
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided stream constructor */
    {
        std::ostringstream ext{};
        sequence_file_output fout{ext, sequence_file_format_fasta{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<sequence_file_format_embl,
                                                                              sequence_file_format_fasta,
                                                                              sequence_file_format_fastq,
                                                                              sequence_file_format_sam>>));// changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided stream temporary constructor */
    {
        sequence_file_output fout{std::ostringstream{}, sequence_file_format_fasta{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<sequence_file_format_embl,
                                                                              sequence_file_format_fasta,
                                                                              sequence_file_format_fastq,
                                                                              sequence_file_format_sam>>));// changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided stream constructor + custom fields + different stream_char_type */
    {
        std::wostringstream ext{};
        sequence_file_output fout{ext, sequence_file_format_fasta{}, fields<field::SEQ>{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));                   // changed
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<sequence_file_format_fasta>>));// changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   wchar_t>));                              // changed
    }

    /* guided stream temporary constructor + custom fields + different stream_char_type */
    {
        sequence_file_output fout{std::wostringstream{}, sequence_file_format_fasta{}, fields<field::SEQ>{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));                   // changed
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<sequence_file_format_fasta>>));// changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   wchar_t>));                              // changed
    }
}

// ----------------------------------------------------------------------------
// *impl
// ----------------------------------------------------------------------------

template <typename fn_t>
void row_wise_impl(fn_t fn)
{
    sequence_file_output fout{std::ostringstream{}, sequence_file_format_fasta{}};
    fout.options.fasta_letters_per_line = 0;

    for (size_t i = 0; i < 3; ++i)
        fn(fout, i);

    fout.get_stream().flush();
    EXPECT_EQ(reinterpret_cast<std::ostringstream&>(fout.get_stream()).str(), output_comp);
}

template <typename source_t>
void assign_impl(source_t && source)
{
    sequence_file_output fout{std::ostringstream{}, sequence_file_format_fasta{}};
    fout.options.fasta_letters_per_line = 0;

    fout = source;

    fout.get_stream().flush();
    EXPECT_EQ(reinterpret_cast<std::ostringstream&>(fout.get_stream()).str(), output_comp);
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

/* Here the record contains a different field composite than the file. The record knows about the
 * association of values and fields, so it does not need to be guessed from the file.
 */
TEST(row, different_fields_in_record_and_file)
{
    std::vector<phred42> qual;
    qual.resize(seqs[1].size());

    record<type_list<std::vector<phred42>, std::string, dna5_vector>,
           fields<field::QUAL, field::ID, field::SEQ>> rec{qual, ids[1], seqs[1]};

    sequence_file_output fout{std::ostringstream{}, sequence_file_format_fasta{}, fields<field::SEQ, field::ID>{}};
    fout.push_back(rec);
    fout.get_stream().flush();

    std::string const expected_out
    {
        "> Test2\n"
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGG\n"
        "CTGNAGGCTGN\n"
    };
    EXPECT_EQ(reinterpret_cast<std::ostringstream&>(fout.get_stream()).str(), expected_out);
}

TEST(row, writing_seq_qual)
{
    sequence_file_output fout{std::ostringstream{}, sequence_file_format_fasta{}, fields<field::ID, field::SEQ_QUAL>()};
    fout.options.fasta_letters_per_line = 0;

    for (size_t i = 0; i < 3; ++i)
    {
        std::vector<qualified<dna5, phred42>> seq_qual;
        seq_qual.resize(quals[i].size());
        std::copy(seqs[i].begin(), seqs[i].end(), seq_qual.begin());
        std::copy(quals[i].begin(), quals[i].end(), seq_qual.begin());

        fout.emplace_back(ids[i], seq_qual);
    }

    EXPECT_EQ(reinterpret_cast<std::ostringstream&>(fout.get_stream()).str(), output_comp);
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

// ----------------------------------------------------------------------------
// columns
// ----------------------------------------------------------------------------

TEST(columns, assign_record_of_columns)
{
    record<type_list<std::vector<dna5_vector>, std::vector<std::string>>, fields<field::SEQ, field::ID> > columns
    {
        seqs,
        ids
    };

    assign_impl(columns);
}

TEST(columns, assign_tuple_of_columns)
{
    assign_impl(std::tie(seqs, ids));
}

TEST(columns, writing_seq_qual)
{
    sequence_file_output fout{std::ostringstream{}, sequence_file_format_fasta{}, fields<field::ID, field::SEQ_QUAL>()};
    fout.options.fasta_letters_per_line = 0;

    std::vector<std::vector<qualified<dna5, phred42>>> seq_quals{3};
    for (size_t i = 0; i < 3; ++i)
    {
        seq_quals[i].resize(quals[i].size());
        std::copy(seqs[i].begin(), seqs[i].end(), seq_quals[i].begin());
        std::copy(quals[i].begin(), quals[i].end(), seq_quals[i].begin());
    }

    fout = std::tie(ids, seq_quals);

    fout.get_stream().flush();
    EXPECT_EQ(reinterpret_cast<std::ostringstream&>(fout.get_stream()).str(), output_comp);
}

// ----------------------------------------------------------------------------
// compression
// ----------------------------------------------------------------------------

void compression_by_filename_impl(test::tmp_filename & filename, std::string_view const expected)
{
    {
        sequence_file_output fout{filename.get_path()};
        fout.options.fasta_letters_per_line = 0;

        for (size_t i = 0; i < 3; ++i)
        {
            record<type_list<dna5_vector, std::string>, fields<field::SEQ, field::ID>> r{seqs[i], ids[i]};

            fout.push_back(r);
        }

    }

    std::string buffer;

    {
        std::ifstream fi{filename.get_path(), std::ios::binary};

        buffer = std::string{std::istreambuf_iterator<char>{fi}, std::istreambuf_iterator<char>{}};
    }

    EXPECT_EQ(buffer, expected);
}

template <typename comp_stream_t>
void compression_by_stream_impl(comp_stream_t & stream)
{
    sequence_file_output fout{stream, sequence_file_format_fasta{}};
    fout.options.fasta_letters_per_line = 0;

    for (size_t i = 0; i < 3; ++i)
    {
        record<type_list<dna5_vector, std::string>, fields<field::SEQ, field::ID>> r{seqs[i], ids[i]};

        fout.push_back(r);
    }
}

#ifdef SEQAN3_HAS_ZLIB
std::string expected_gz
{
    '\x1F','\x8B','\x08','\x00','\x00','\x00','\x00','\x00','\x00','\x03','\xB3','\x53','\x08','\x71','\x0D','\x0E',
    '\x51','\x30','\xE4','\x72','\x74','\x76','\x0F','\xE1','\xB2','\x53','\x08','\x49','\x2D','\x2E','\x31','\xE2',
    '\x72','\x74','\x77','\x77','\x0E','\x71','\xF7','\xA3','\x05','\x05','\xB5','\xC3','\x98','\xCB','\xDD','\xDD',
    '\xD1','\x3D','\xC4','\x31','\xC4','\xD1','\x31','\x04','\x15','\x72','\x01','\x00','\x27','\xAD','\xB4','\xE9',
    '\x93','\x00','\x00','\x00'
};

TEST(compression, by_filename_gz)
{
    test::tmp_filename filename{"sequence_file_output_test.fasta.gz"};

    compression_by_filename_impl(filename, expected_gz);
}

TEST(compression, by_stream_gz)
{
    std::ostringstream out;

    {
        contrib::gz_ostream compout{out};
        compression_by_stream_impl(compout);
    }

    EXPECT_EQ(out.str(), expected_gz);
}
#endif

#ifdef SEQAN3_HAS_BZIP2
std::string expected_bz2
{
    '\x42','\x5A','\x68','\x39','\x31','\x41','\x59','\x26','\x53','\x59','\xB4','\x68','\xEA','\xE3','\x00','\x00',
    '\x06','\xDF','\x80','\x00','\x10','\x40','\x00','\x38','\x01','\x2A','\x81','\x0C','\x00','\x02','\x00','\x0C',
    '\x00','\x20','\x00','\x50','\xA6','\x00','\x09','\xA0','\x8A','\x10','\x9A','\x32','\x34','\xD9','\xAB','\x5F',
    '\x16','\xE9','\xEB','\x86','\x5B','\x46','\x41','\x8D','\xD0','\x1E','\x12','\x8C','\xC0','\xB5','\x48','\xD2',
    '\x3A','\x9B','\x23','\xB9','\x9F','\x64','\x98','\x1E','\xEE','\x8C','\x18','\x3E','\x38','\x7E','\x2E','\xE4',
    '\x8A','\x70','\xA1','\x21','\x68','\xD1','\xD5','\xC6'
};

TEST(compression, by_filename_bz2)
{
    test::tmp_filename filename{"sequence_file_output_test.fasta.bz2"};


    compression_by_filename_impl(filename, expected_bz2);
}

TEST(compression, by_stream_bz2)
{
    std::ostringstream out;

    {
        contrib::bz2_ostream compout{out};
        compression_by_stream_impl(compout);
    }

    EXPECT_EQ(out.str(), expected_bz2);
}
#endif
