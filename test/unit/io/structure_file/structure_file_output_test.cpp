// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <iterator>
#include <fstream>
#include <sstream>
#include <vector>

#include <gtest/gtest.h>

#include <range/v3/view/zip.hpp>
#include <range/v3/view/filter.hpp>
#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/map.hpp>

#include <seqan3/io/structure_file/output.hpp>
#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <seqan3/test/tmp_filename.hpp>

using namespace seqan3;

TEST(general, concepts)
{
    using it_t = typename structure_file_out<>::iterator;
    using sen_t = typename structure_file_out<>::sentinel;

    EXPECT_TRUE((std::OutputIterator<it_t, std::tuple<std::string, std::string>>));
    EXPECT_TRUE((std::Sentinel<sen_t, it_t>));
}

TEST(structure_file_output_class, concepts)
{
    using t = structure_file_out<>;
    EXPECT_TRUE((std::ranges::OutputRange<t, std::tuple<std::string, std::string>>));

    using ct = structure_file_out<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::OutputRange<ct, std::tuple<std::string, std::string>>));
}

TEST(structure_file_output_class, construct_by_filename)
{
    /* just the filename */
    {
        test::tmp_filename filename{"structure_file_output_constructor.dbn"};
        EXPECT_NO_THROW(structure_file_out<>{filename.get_path()});
    }

    /* wrong extension */
    {
        test::tmp_filename filename{"structure_file_output_constructor.xyz"};
        EXPECT_THROW(structure_file_out<>{filename.get_path()}, unhandled_extension_error);
    }

    /* unknown file */
    {
        test::tmp_filename filename{"I/do/not/exist.dbn"};
        EXPECT_THROW(structure_file_out<>{filename.get_path()}, file_open_error);
    }

    /* non-existent file*/
    {
        EXPECT_THROW(structure_file_out<>{"/dev/nonexistant/foobarOOO"}, file_open_error);
    }

    /* filename + fields */
    {
        test::tmp_filename filename{"structure_file_output_constructor.dbn"};
        EXPECT_NO_THROW((structure_file_out<fields<field::SEQ>,
                                            type_list<structure_file_format_vienna>>
                                            {filename.get_path(), fields<field::SEQ>{}}));
    }
}

TEST(structure_file_output_class, construct_from_stream)
{
    /* stream + format_tag */
    EXPECT_NO_THROW((structure_file_out<fields<field::SEQ, field::ID, field::STRUCTURE>,
                                        type_list<structure_file_format_vienna>>
                                        {std::ostringstream{}, structure_file_format_vienna{}}));

    /* stream + format_tag + fields */
    EXPECT_NO_THROW((structure_file_out<fields<field::SEQ, field::ID, field::STRUCTURE>,
                                        type_list<structure_file_format_vienna>>
                     {std::ostringstream{}, structure_file_format_vienna{},
                      fields<field::SEQ, field::ID, field::STRUCTURE>{}}));
}

TEST(structure_file_output_class, default_template_args_and_deduction_guides)
{
    using comp1 = fields<field::SEQ, field::ID, field::STRUCTURE>;
    using comp2 = type_list<structure_file_format_vienna>;
    using comp3 = char;

    /* default template args */
    {
        using t = structure_file_out<>;
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided filename constructor */
    {
        test::tmp_filename filename{"structure_file_output_constructor.dbn"};
        structure_file_out fout{filename.get_path()};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided filename constructor + custom fields */
    {
        test::tmp_filename filename{"structure_file_output_constructor.dbn"};
        structure_file_out fout{filename.get_path(), fields<field::SEQ>{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided stream constructor */
    {
        std::ostringstream ext{};
        structure_file_out fout{ext, structure_file_format_vienna{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<structure_file_format_vienna>>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided stream temporary constructor */
    {
        structure_file_out fout{std::ostringstream{}, structure_file_format_vienna{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<structure_file_format_vienna>>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided stream constructor + custom fields + different stream_char_type */
    {
        std::wostringstream ext{};
        structure_file_out fout{ext, structure_file_format_vienna{}, fields<field::SEQ>{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<structure_file_format_vienna>>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   wchar_t>));
    }

    /* guided stream temporary constructor + custom fields + different stream_char_type */
    {
        structure_file_out fout{std::wostringstream{}, structure_file_format_vienna{}, fields<field::SEQ>{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<structure_file_format_vienna>>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   wchar_t>));
    }
}

struct structure_file_output_write : public ::testing::Test
{
    size_t const num_records = 2ul;

    std::vector<rna5_vector> const seqs
    {
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA"_rna5,
        "UUGGAGUACACAACCUGUACACUCUUUC"_rna5
    };

    std::vector<std::string> const ids
    {
        "S.cerevisiae_tRNA-PHE M10740/1-73",
        "example"
    };

    std::vector<double> const energies
    {
        -17.5, -3.71
    };

    std::vector<std::vector<wuss51>> const structures
    {
        "(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))."_wuss51,
        "..(((((..(((...)))..)))))..."_wuss51
    };

    std::string const output_comp
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).\n"
        "> example\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))...\n"
    };
};

struct structure_file_output_row : public structure_file_output_write
{
    template <typename fn_t>
    void row_wise_impl(fn_t fn)
    {
        structure_file_out fout{std::ostringstream{}, structure_file_format_vienna{}};

        for (size_t idx = 0ul; idx < num_records; ++idx)
            fn(fout, idx);

        fout.get_stream().flush();
        EXPECT_EQ(reinterpret_cast<std::ostringstream&>(fout.get_stream()).str(), output_comp);
    }
};

TEST_F(structure_file_output_row, assign_to_iterator)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        record<type_list<rna5_vector, std::string, std::vector<wuss51>>,
               fields<field::SEQ, field::ID, field::STRUCTURE>> r{seqs[i], ids[i], structures[i]};
        begin(file) = r;
    });
}

TEST_F(structure_file_output_row, push_back_record)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        record<type_list<rna5_vector, std::string, std::vector<wuss51>>,
               fields<field::SEQ, field::ID, field::STRUCTURE>> r{seqs[i], ids[i], structures[i]};
        file.push_back(r);
    });
}

TEST_F(structure_file_output_row, push_back_record_rvalue)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        record<type_list<rna5_vector, std::string, std::vector<wuss51>>,
               fields<field::SEQ, field::ID, field::STRUCTURE>> r{seqs[i], ids[i], structures[i]};
        file.push_back(std::move(r));
    });
}

TEST_F(structure_file_output_row, push_back_record_const)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        record<type_list<rna5_vector, std::string, std::vector<wuss51>>,
               fields<field::SEQ, field::ID, field::STRUCTURE>> const r{seqs[i], ids[i], structures[i]};
        file.push_back(r);
    });
}

TEST_F(structure_file_output_row, push_back_record_const_element)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        record<type_list<rna5_vector const, std::string const, std::vector<wuss51> const>,
               fields<field::SEQ, field::ID, field::STRUCTURE>> const r{seqs[i], ids[i], structures[i]};
        file.push_back(r);
    });
}

TEST_F(structure_file_output_row, push_back_tuple)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        std::tuple<rna5_vector, std::string, std::vector<wuss51>> t{seqs[i], ids[i], structures[i]};
        file.push_back(t);
    });
}

TEST_F(structure_file_output_row, push_back_tuple_rvalue)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        std::tuple<rna5_vector, std::string, std::vector<wuss51>> t{seqs[i], ids[i], structures[i]};
        file.push_back(std::move(t));
    });
}

TEST_F(structure_file_output_row, push_back_tuple_const)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        std::tuple<rna5_vector, std::string, std::vector<wuss51>> const t{seqs[i], ids[i], structures[i]};
        file.push_back(t);
    });
}

TEST_F(structure_file_output_row, push_back_tuple_const_element)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        std::tuple<rna5_vector const, std::string const, std::vector<wuss51> const> t{seqs[i], ids[i], structures[i]};
        file.push_back(t);
    });
}

TEST_F(structure_file_output_row, emplace_back)
{
    row_wise_impl([&] (auto & file, size_t i)
    {
        file.emplace_back(seqs[i], ids[i], structures[i]);
    });
}

struct structure_file_output_rows : public structure_file_output_write
{
    template<typename source_t>
    void assign_impl(source_t && source)
    {
        structure_file_out fout{std::ostringstream{}, structure_file_format_vienna{}};
        fout = source;
        fout.get_stream().flush();
        EXPECT_EQ(reinterpret_cast<std::ostringstream&>(fout.get_stream()).str(), output_comp);
    }
};

TEST_F(structure_file_output_rows, assign_range_of_records)
{
    std::vector<record<type_list<rna5_vector, std::string, std::vector<wuss51>>,
                fields<field::SEQ, field::ID, field::STRUCTURE>>> range;

    for (size_t idx = 0ul; idx < num_records; ++idx)
        range.emplace_back(seqs[idx], ids[idx], structures[idx]);

    assign_impl(range);
}

TEST_F(structure_file_output_rows, assign_range_of_records_const)
{
    std::vector<record<type_list<rna5_vector, std::string, std::vector<wuss51>>,
                fields<field::SEQ, field::ID, field::STRUCTURE>>> range;

    for (size_t idx = 0ul; idx < num_records; ++idx)
        range.emplace_back(seqs[idx], ids[idx], structures[idx]);

    assign_impl(std::as_const(range));
}

TEST_F(structure_file_output_rows, assign_range_of_tuples)
{
    std::vector<std::tuple<rna5_vector, std::string, std::vector<wuss51>>> range;

    for (size_t idx = 0ul; idx < num_records; ++idx)
        range.emplace_back(seqs[idx], ids[idx], structures[idx]);

    assign_impl(range);
}

TEST_F(structure_file_output_rows, assign_structure_file_in)
{
    std::string const inp // differs from output above by formatting
    {
        "> S.cerevisiae_tRNA-PHE M10740/1-73\n"
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
        "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).\n"
        "> example\n"
        "UUGGAGUACACAACCUGUACACUCUUUC\n"
        "..(((((..(((...)))..)))))...\n"
    };

    structure_file_in fin{std::istringstream{inp}, structure_file_format_vienna{},
                          fields<field::SEQ, field::ID, field::STRUCTURE>{}};
    assign_impl(fin);
}

TEST_F(structure_file_output_rows, assign_structure_file_pipes)
{
    // valid without assignment?
    structure_file_in{std::istringstream{output_comp}, structure_file_format_vienna{}}
              | structure_file_out{std::ostringstream{}, structure_file_format_vienna{}};

    // valid with assignment and check contents
    auto fout = structure_file_in{std::istringstream{output_comp}, structure_file_format_vienna{}}
              | structure_file_out{std::ostringstream{}, structure_file_format_vienna{}};

    fout.get_stream().flush();
    EXPECT_EQ(reinterpret_cast<std::ostringstream&>(fout.get_stream()).str(), output_comp);
}

struct structure_file_output_columns : public structure_file_output_rows{};

TEST_F(structure_file_output_columns, assign_record_of_columns)
{
    record<type_list<std::vector<rna5_vector>, std::vector<std::string>, std::vector<std::vector<wuss51>>>,
           fields<field::SEQ, field::ID, field::STRUCTURE>> columns
    {
        seqs,
        ids,
        structures
    };

    assign_impl(columns);
}

TEST_F(structure_file_output_columns, assign_tuple_of_columns)
{
    assign_impl(std::tie(seqs, ids, structures));
}

// ----------------------------------------------------------------------------
// compression
// ----------------------------------------------------------------------------

#if defined(SEQAN3_HAS_ZLIB) || defined(SEQAN3_HAS_BZIP2)
struct structure_file_output_compression : public structure_file_output_write
{
    void compression_by_filename_impl(test::tmp_filename & filename, std::string_view const expected)
    {
        {
            structure_file_out
            fout{filename.get_path()};

            for (size_t idx = 0ul; idx < num_records; ++idx)
            {
                record<type_list<rna5_vector, std::string, std::vector<wuss51>>,
                       fields<field::SEQ, field::ID, field::STRUCTURE>> rec{seqs[idx], ids[idx], structures[idx]};
                fout.push_back(rec);
            }
        }
        std::string buffer;
        {
            std::ifstream fi{filename.get_path(), std::ios::binary};
            buffer = std::string{std::istreambuf_iterator<char>{fi}, std::istreambuf_iterator<char>{}};
        }
        EXPECT_EQ(buffer, expected);
    }

    template<typename comp_stream_t>
    void compression_by_stream_impl(comp_stream_t & stream)
    {
        structure_file_out fout{stream, structure_file_format_vienna{}};
        for (size_t idx = 0ul; idx < num_records; ++idx)
        {
            record<type_list<rna5_vector, std::string, std::vector<wuss51>>,
                   fields<field::SEQ, field::ID, field::STRUCTURE>> rec{seqs[idx], ids[idx], structures[idx]};
            fout.push_back(rec);
        }
    }
};
#endif

#ifdef SEQAN3_HAS_ZLIB
std::string expected_gz
{
    '\x1F','\x8B','\x08','\x00','\x00','\x00','\x00','\x00','\x00','\x03','\x55','\x8E','\xC1','\x0A','\xC2','\x40',
    '\x0C','\x44','\xEF','\xF9','\x8A','\x3D','\x76','\x0F','\x5D','\x5B','\x14','\x7A','\x2B','\x84','\x20','\xF1',
    '\xA2','\x88','\x92','\xB3','\x14','\xD9','\x43','\x41','\x41','\xB4','\x14','\x3F','\xDF','\x64','\x23','\x52',
    '\x27','\xB0','\x64','\x1E','\x61','\x66','\xFB','\x70','\x4E','\xD7','\xFC','\xCC','\xF3','\xF8','\x1A','\x87',
    '\x7C','\x99','\x4E','\x07','\xAC','\x8F','\xBB','\x6D','\xD8','\xB7','\x4D','\xB7','\x69','\x56','\x6D','\xDD',
    '\xAD','\x81','\x89','\x19','\x45','\x04','\x99','\x84','\x90','\x45','\x58','\xBD','\x0E','\x31','\xA9','\x45',
    '\x12','\x46','\x2C','\x07','\x86','\x59','\x48','\x81','\x8E','\x90','\x32','\x3D','\xB0','\x13','\x34','\x47',
    '\x08','\x95','\x2B','\x25','\x7F','\x5D','\x51','\xF5','\x07','\x9C','\x98','\xAA','\x05','\x8E','\x0B','\x25',
    '\xE8','\x43','\x7E','\x0F','\xF7','\xC7','\x2D','\x83','\xD7','\x0A','\x5A','\x13','\x96','\x6E','\x5B','\x85',
    '\xF4','\x3F','\x04','\xBF','\x08','\xCF','\x29','\xB9','\xF1','\x1B','\x0F','\x1F','\xA0','\x5A','\xBE','\x54',
    '\xFC','\x00','\x00','\x00'
};

TEST_F(structure_file_output_compression, by_filename_gz)
{
    test::tmp_filename filename{"structure_file_output_test.dbn.gz"};
    compression_by_filename_impl(filename, expected_gz);
}

TEST_F(structure_file_output_compression, by_stream_gz)
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
    '\x42','\x5A','\x68','\x39','\x31','\x41','\x59','\x26','\x53','\x59','\xC7','\x0B','\xB5','\x7F','\x00','\x00',
    '\x36','\x5F','\x80','\x6E','\x50','\x40','\x63','\xEC','\x81','\x2A','\xC3','\x5A','\x00','\xAA','\x26','\x5D',
    '\x40','\x30','\x00','\xB8','\x84','\x53','\xC5','\x00','\x68','\x00','\x03','\x40','\x34','\x69','\xEA','\x18',
    '\x01','\x93','\x4D','\x06','\x43','\x04','\x34','\xC4','\x68','\xC0','\x94','\xD3','\x52','\xA7','\xEA','\x9B',
    '\x14','\xF2','\x69','\x1E','\xA7','\xA9','\xEA','\x68','\xC0','\x23','\x4D','\x35','\x85','\x85','\xCA','\x54',
    '\xA4','\x4F','\xB6','\x4C','\xD9','\xCB','\x3C','\xCD','\x51','\x11','\xE5','\x16','\xEB','\x96','\x5A','\x11',
    '\x7E','\x14','\xC1','\x50','\xCB','\x07','\x06','\x2B','\x15','\x01','\x5B','\x6E','\xD5','\x48','\x26','\xEA',
    '\xCA','\x37','\x7B','\xE7','\xE9','\x9E','\xDD','\x0D','\x2B','\x79','\xF1','\xF4','\xB6','\x8B','\x78','\xB2',
    '\x4D','\x0A','\x53','\x43','\x4D','\x0D','\x48','\xD0','\x98','\xDC','\xC4','\xC4','\x8C','\x7F','\x69','\x94',
    '\x48','\xA2','\x99','\x15','\x53','\xA1','\x44','\xC1','\x31','\x02','\x5A','\xF5','\x91','\xA7','\x00','\x40',
    '\x10','\xC2','\x66','\x06','\x02','\xE0','\x81','\x10','\x09','\x94','\x46','\x6E','\x8E','\xBD','\x26','\x2C',
    '\xED','\x8D','\x97','\xE4','\x47','\xD1','\x4A','\x42','\x0F','\xC5','\xDC','\x91','\x4E','\x14','\x24','\x31',
    '\xC2','\xED','\x5F','\xC0'
};

TEST_F(structure_file_output_compression, by_filename_bz2)
{
    test::tmp_filename filename{"structure_file_output_test.dbn.bz2"};
    compression_by_filename_impl(filename, expected_bz2);
}

TEST_F(structure_file_output_compression, by_stream_bz2)
{
    std::ostringstream out;
    {
        contrib::bz2_ostream compout{out};
        compression_by_stream_impl(compout);
    }
    EXPECT_EQ(out.str(), expected_bz2);
}
#endif
