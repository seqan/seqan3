// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <fstream>
#include <iterator>
#include <ranges>
#include <sstream>
#include <vector>

#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/io/structure_file/output.hpp>
#include <seqan3/test/tmp_directory.hpp>
#include <seqan3/test/zlib_skip.hpp>

using seqan3::operator""_rna5;
using seqan3::operator""_wuss51;

TEST(general, concepts)
{
    using it_t = typename seqan3::structure_file_output<>::iterator;
    using sen_t = typename seqan3::structure_file_output<>::sentinel;

    EXPECT_TRUE((std::output_iterator<it_t, std::tuple<std::string, std::string>>));
    EXPECT_TRUE((std::sentinel_for<sen_t, it_t>));
}

TEST(structure_file_output_class, concepts)
{
    using t = seqan3::structure_file_output<>;
    EXPECT_TRUE((std::ranges::output_range<t, std::tuple<std::string, std::string>>));

    using ct = seqan3::structure_file_output<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::output_range<ct, std::tuple<std::string, std::string>>));
}

TEST(structure_file_output_class, construct_by_filename)
{
    /* just the filename */
    {
        seqan3::test::tmp_directory tmp;
        auto filename = tmp.path() / "structure_file_output_constructor.dbn";
        EXPECT_NO_THROW(seqan3::structure_file_output<>{filename});
    }

    /* wrong extension */
    {
        seqan3::test::tmp_directory tmp;
        auto filename = tmp.path() / "structure_file_output_constructor.xyz";
        EXPECT_THROW(seqan3::structure_file_output<>{filename}, seqan3::unhandled_extension_error);
    }

    /* unknown file */
    {
        seqan3::test::tmp_directory tmp;
        auto filename = tmp.path() / "I/do/not/exist.dbn";
        EXPECT_THROW(seqan3::structure_file_output<>{filename}, seqan3::file_open_error);
    }

    /* non-existent file*/
    {
        EXPECT_THROW(seqan3::structure_file_output<>{"/dev/nonexistant/foobarOOO"}, seqan3::file_open_error);
    }

    /* filename + fields */
    {
        seqan3::test::tmp_directory tmp;
        auto filename = tmp.path() / "structure_file_output_constructor.dbn";
        EXPECT_NO_THROW((
            seqan3::structure_file_output<seqan3::fields<seqan3::field::seq>, seqan3::type_list<seqan3::format_vienna>>{
                filename,
                seqan3::fields<seqan3::field::seq>{}}));
    }
}

TEST(structure_file_output_class, construct_from_stream)
{
    using fields_seq_id_structure = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>;

    /* stream + format_tag */
    EXPECT_NO_THROW((seqan3::structure_file_output<fields_seq_id_structure, seqan3::type_list<seqan3::format_vienna>>{
        std::ostringstream{},
        seqan3::format_vienna{}}));

    /* stream + format_tag + fields */
    EXPECT_NO_THROW((seqan3::structure_file_output<fields_seq_id_structure, seqan3::type_list<seqan3::format_vienna>>{
        std::ostringstream{},
        seqan3::format_vienna{},
        fields_seq_id_structure{}}));
}

TEST(structure_file_output_class, default_template_args_and_deduction_guides)
{
    using comp1 = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>;
    using comp2 = seqan3::type_list<seqan3::format_vienna>;
    using comp3 = char;

    /* default template args */
    {
        using t = seqan3::structure_file_output<>;
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats, comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
    }

    /* guided filename constructor */
    {
        seqan3::test::tmp_directory tmp;
        auto filename = tmp.path() / "structure_file_output_constructor.dbn";
        seqan3::structure_file_output fout{filename};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats, comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
    }

    /* guided filename constructor + custom fields */
    {
        seqan3::test::tmp_directory tmp;
        auto filename = tmp.path() / "structure_file_output_constructor.dbn";
        seqan3::structure_file_output fout{filename, seqan3::fields<seqan3::field::seq>{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, seqan3::fields<seqan3::field::seq>>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats, comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
    }

    /* guided stream constructor */
    {
        std::ostringstream ext{};
        seqan3::structure_file_output fout{ext, seqan3::format_vienna{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats, seqan3::type_list<seqan3::format_vienna>>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
    }

    /* guided stream temporary constructor */
    {
        seqan3::structure_file_output fout{std::ostringstream{}, seqan3::format_vienna{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats, seqan3::type_list<seqan3::format_vienna>>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
    }
}

struct structure_file_output_write : public ::testing::Test
{
    size_t const num_records = 2ul;

    std::vector<seqan3::rna5_vector> const seqs{
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA"_rna5,
        "UUGGAGUACACAACCUGUACACUCUUUC"_rna5};

    std::vector<std::string> const ids{"S.cerevisiae_tRNA-PHE M10740/1-73", "example"};

    std::vector<double> const energies{-17.5, -3.71};

    std::vector<std::vector<seqan3::wuss51>> const structures{
        "(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))."_wuss51,
        "..(((((..(((...)))..)))))..."_wuss51};

    std::string const output_comp{"> S.cerevisiae_tRNA-PHE M10740/1-73\n"
                                  "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
                                  "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).\n"
                                  "> example\n"
                                  "UUGGAGUACACAACCUGUACACUCUUUC\n"
                                  "..(((((..(((...)))..)))))...\n"};
};

struct structure_file_output_row : public structure_file_output_write
{
    template <typename fn_t>
    void row_wise_impl(fn_t fn)
    {
        seqan3::structure_file_output fout{std::ostringstream{}, seqan3::format_vienna{}};

        for (size_t idx = 0ul; idx < num_records; ++idx)
            fn(fout, idx);

        fout.get_stream().flush();
        EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), output_comp);
    }
};

TEST_F(structure_file_output_row, assign_to_iterator)
{
    row_wise_impl(
        [&](auto & file, size_t i)
        {
            seqan3::record<seqan3::type_list<seqan3::rna5_vector, std::string, std::vector<seqan3::wuss51>>,
                           seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>>
                r{seqs[i], ids[i], structures[i]};
            std::ranges::begin(file) = r;
        });
}

TEST_F(structure_file_output_row, push_back_record)
{
    row_wise_impl(
        [&](auto & file, size_t i)
        {
            seqan3::record<seqan3::type_list<seqan3::rna5_vector, std::string, std::vector<seqan3::wuss51>>,
                           seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>>
                r{seqs[i], ids[i], structures[i]};
            file.push_back(r);
        });
}

TEST_F(structure_file_output_row, push_back_record_rvalue)
{
    row_wise_impl(
        [&](auto & file, size_t i)
        {
            seqan3::record<seqan3::type_list<seqan3::rna5_vector, std::string, std::vector<seqan3::wuss51>>,
                           seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>>
                r{seqs[i], ids[i], structures[i]};
            file.push_back(std::move(r));
        });
}

TEST_F(structure_file_output_row, push_back_record_const)
{
    row_wise_impl(
        [&](auto & file, size_t i)
        {
            seqan3::record<seqan3::type_list<seqan3::rna5_vector, std::string, std::vector<seqan3::wuss51>>,
                           seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>> const r{
                seqs[i],
                ids[i],
                structures[i]};
            file.push_back(r);
        });
}

TEST_F(structure_file_output_row, push_back_record_const_element)
{
    row_wise_impl(
        [&](auto & file, size_t i)
        {
            seqan3::record<
                seqan3::type_list<seqan3::rna5_vector const, std::string const, std::vector<seqan3::wuss51> const>,
                seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>> const r{seqs[i],
                                                                                                         ids[i],
                                                                                                         structures[i]};
            file.push_back(r);
        });
}

TEST_F(structure_file_output_row, push_back_tuple)
{
    row_wise_impl(
        [&](auto & file, size_t i)
        {
            std::tuple<seqan3::rna5_vector, std::string, std::vector<seqan3::wuss51>> t{seqs[i], ids[i], structures[i]};
            file.push_back(t);
        });
}

TEST_F(structure_file_output_row, push_back_tuple_rvalue)
{
    row_wise_impl(
        [&](auto & file, size_t i)
        {
            std::tuple<seqan3::rna5_vector, std::string, std::vector<seqan3::wuss51>> t{seqs[i], ids[i], structures[i]};
            file.push_back(std::move(t));
        });
}

TEST_F(structure_file_output_row, push_back_tuple_const)
{
    row_wise_impl(
        [&](auto & file, size_t i)
        {
            std::tuple<seqan3::rna5_vector, std::string, std::vector<seqan3::wuss51>> const t{seqs[i],
                                                                                              ids[i],
                                                                                              structures[i]};
            file.push_back(t);
        });
}

TEST_F(structure_file_output_row, push_back_tuple_const_element)
{
    row_wise_impl(
        [&](auto & file, size_t i)
        {
            std::tuple<seqan3::rna5_vector const, std::string const, std::vector<seqan3::wuss51> const> t{
                seqs[i],
                ids[i],
                structures[i]};
            file.push_back(t);
        });
}

TEST_F(structure_file_output_row, emplace_back)
{
    row_wise_impl(
        [&](auto & file, size_t i)
        {
            file.emplace_back(seqs[i], ids[i], structures[i]);
        });
}

struct structure_file_output_rows : public structure_file_output_write
{
    template <typename source_t>
    void assign_impl(source_t && source)
    {
        seqan3::structure_file_output fout{std::ostringstream{}, seqan3::format_vienna{}};
        fout = source;
        fout.get_stream().flush();
        EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), output_comp);
    }
};

TEST_F(structure_file_output_rows, assign_range_of_records)
{
    std::vector<seqan3::record<seqan3::type_list<seqan3::rna5_vector, std::string, std::vector<seqan3::wuss51>>,
                               seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>>>
        range;

    for (size_t idx = 0ul; idx < num_records; ++idx)
        range.emplace_back(seqs[idx], ids[idx], structures[idx]);

    assign_impl(range);
}

TEST_F(structure_file_output_rows, assign_range_of_records_const)
{
    std::vector<seqan3::record<seqan3::type_list<seqan3::rna5_vector, std::string, std::vector<seqan3::wuss51>>,
                               seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>>>
        range;

    for (size_t idx = 0ul; idx < num_records; ++idx)
        range.emplace_back(seqs[idx], ids[idx], structures[idx]);

    assign_impl(std::as_const(range));
}

TEST_F(structure_file_output_rows, assign_range_of_tuples)
{
    std::vector<std::tuple<seqan3::rna5_vector, std::string, std::vector<seqan3::wuss51>>> range;

    for (size_t idx = 0ul; idx < num_records; ++idx)
        range.emplace_back(seqs[idx], ids[idx], structures[idx]);

    assign_impl(range);
}

TEST_F(structure_file_output_rows, assign_structure_file_input)
{
    std::string const inp // differs from output above by formatting
        {"> S.cerevisiae_tRNA-PHE M10740/1-73\n"
         "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
         "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))).\n"
         "> example\n"
         "UUGGAGUACACAACCUGUACACUCUUUC\n"
         "..(((((..(((...)))..)))))...\n"};

    seqan3::structure_file_input fin{std::istringstream{inp},
                                     seqan3::format_vienna{},
                                     seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>{}};
    assign_impl(fin);
}

TEST_F(structure_file_output_rows, assign_structure_file_pipes)
{
    // valid without assignment?
    seqan3::structure_file_input{std::istringstream{output_comp}, seqan3::format_vienna{}}
        | seqan3::structure_file_output{std::ostringstream{}, seqan3::format_vienna{}};

    // valid with assignment and check contents
    auto fout = seqan3::structure_file_input{std::istringstream{output_comp}, seqan3::format_vienna{}}
              | seqan3::structure_file_output{std::ostringstream{}, seqan3::format_vienna{}};

    fout.get_stream().flush();
    EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), output_comp);
}

struct structure_file_output_columns : public structure_file_output_rows
{};

TEST_F(structure_file_output_columns, assign_columns)
{
    assign_impl(seqan3::views::zip(seqs, ids, structures));
}

// ----------------------------------------------------------------------------
// compression
// ----------------------------------------------------------------------------

#if defined(SEQAN3_HAS_ZLIB) || defined(SEQAN3_HAS_BZIP2)
struct structure_file_output_compression : public structure_file_output_write
{
    std::string compression_by_filename_impl(seqan3::test::sandboxed_path const & filename)
    {
        {
            seqan3::structure_file_output fout{filename};

            for (size_t idx = 0ul; idx < num_records; ++idx)
            {
                seqan3::record<seqan3::type_list<seqan3::rna5_vector, std::string, std::vector<seqan3::wuss51>>,
                               seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>>
                    rec{seqs[idx], ids[idx], structures[idx]};
                fout.push_back(rec);
            }
        }
        std::string buffer{};
        {
            std::ifstream fi{filename, std::ios::binary};
            buffer = std::string{std::istreambuf_iterator<char>{fi}, std::istreambuf_iterator<char>{}};
        }
        return buffer;
    }

    template <typename comp_stream_t>
    void compression_by_stream_impl(comp_stream_t & stream)
    {
        seqan3::structure_file_output fout{stream, seqan3::format_vienna{}};
        for (size_t idx = 0ul; idx < num_records; ++idx)
        {
            seqan3::record<seqan3::type_list<seqan3::rna5_vector, std::string, std::vector<seqan3::wuss51>>,
                           seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>>
                rec{seqs[idx], ids[idx], structures[idx]};
            fout.push_back(rec);
        }
    }
};
#endif

#if defined(SEQAN3_HAS_ZLIB)
std::string expected_gz{
    '\x1F', '\x8B', '\x08', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x55', '\x8E', '\xC1', '\x0A',
    '\xC2', '\x40', '\x0C', '\x44', '\xEF', '\xF9', '\x8A', '\x3D', '\x76', '\x0F', '\x5D', '\x5B', '\x14', '\x7A',
    '\x2B', '\x84', '\x20', '\xF1', '\xA2', '\x88', '\x92', '\xB3', '\x14', '\xD9', '\x43', '\x41', '\x41', '\xB4',
    '\x14', '\x3F', '\xDF', '\x64', '\x23', '\x52', '\x27', '\xB0', '\x64', '\x1E', '\x61', '\x66', '\xFB', '\x70',
    '\x4E', '\xD7', '\xFC', '\xCC', '\xF3', '\xF8', '\x1A', '\x87', '\x7C', '\x99', '\x4E', '\x07', '\xAC', '\x8F',
    '\xBB', '\x6D', '\xD8', '\xB7', '\x4D', '\xB7', '\x69', '\x56', '\x6D', '\xDD', '\xAD', '\x81', '\x89', '\x19',
    '\x45', '\x04', '\x99', '\x84', '\x90', '\x45', '\x58', '\xBD', '\x0E', '\x31', '\xA9', '\x45', '\x12', '\x46',
    '\x2C', '\x07', '\x86', '\x59', '\x48', '\x81', '\x8E', '\x90', '\x32', '\x3D', '\xB0', '\x13', '\x34', '\x47',
    '\x08', '\x95', '\x2B', '\x25', '\x7F', '\x5D', '\x51', '\xF5', '\x07', '\x9C', '\x98', '\xAA', '\x05', '\x8E',
    '\x0B', '\x25', '\xE8', '\x43', '\x7E', '\x0F', '\xF7', '\xC7', '\x2D', '\x83', '\xD7', '\x0A', '\x5A', '\x13',
    '\x96', '\x6E', '\x5B', '\x85', '\xF4', '\x3F', '\x04', '\xBF', '\x08', '\xCF', '\x29', '\xB9', '\xF1', '\x1B',
    '\x0F', '\x1F', '\xA0', '\x5A', '\xBE', '\x54', '\xFC', '\x00', '\x00', '\x00'};

std::string expected_bgzf{
    '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x06', '\x00', '\x42', '\x43',
    '\x02', '\x00', '\xAF', '\x00', '\x55', '\x4E', '\xB1', '\x0A', '\xC2', '\x50', '\x0C', '\xDC', '\xF3', '\x15',
    '\x1D', '\x75', '\x68', '\x6C', '\x51', '\xE8', '\x56', '\x08', '\x41', '\xE2', '\xA2', '\x88', '\x92', '\x59',
    '\x8A', '\xBC', '\xA1', '\xA0', '\x50', '\x54', '\xC4', '\xCF', '\xF7', '\xF2', '\x9E', '\x88', '\x5E', '\x86',
    '\x7B', '\xB9', '\x5C', '\x72', '\xAF', '\xAF', '\x8E', '\x7C', '\x4E', '\xB7', '\xF4', '\x1C', '\xEF', '\xE3',
    '\x90', '\x4E', '\x8F', '\xC3', '\x4E', '\xEA', '\xFD', '\x66', '\x5D', '\x6D', '\xDB', '\xA6', '\x5B', '\x35',
    '\x8B', '\xB6', '\xEE', '\x96', '\x64', '\x6A', '\x26', '\xEE', '\x2E', '\xA6', '\xAE', '\x62', '\xEE', '\x86',
    '\x1E', '\xA5', '\xA6', '\x68', '\x45', '\xDD', '\x04', '\x04', '\x43', '\xC8', '\xE6', '\x0A', '\x01', '\xE5',
    '\x0A', '\x0D', '\x86', '\xB0', '\x60', '\x08', '\xB3', '\xD0', '\xAC', '\x80', '\x39', '\x98', '\x3F', '\x98',
    '\x03', '\x7F', '\x02', '\x67', '\x25', '\xA6', '\xD9', '\xFE', '\x63', '\x8B', '\x41', '\x80', '\xA9', '\xAF',
    '\xD2', '\x6B', '\xB8', '\x4E', '\x97', '\x44', '\x25', '\xD6', '\x91', '\xA3', '\x22', '\x39', '\x3B', '\x9E',
    '\xAE', '\xF8', '\x8F', '\xD2', '\xF7', '\x44', '\xC9', '\x8B', '\xD5', '\x7C', '\x1D', '\xC4', '\xF4', '\x06',
    '\xA0', '\x5A', '\xBE', '\x54', '\xFC', '\x00', '\x00', '\x00', '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00',
    '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43', '\x02', '\x00', '\x1B', '\x00', '\x03', '\x00',
    '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00'};

TEST_F(structure_file_output_compression, by_filename_gz)
{
    seqan3::test::tmp_directory tmp;
    auto filename = tmp.path() / "structure_file_output_test.dbn.gz";

    std::string buffer = compression_by_filename_impl(filename);
    buffer[9] = '\x00'; // zero out OS byte
    SEQAN3_TEST_GTEST_SKIP_ZLIB_DEFLATE;
    EXPECT_EQ(buffer, expected_gz);
}

TEST_F(structure_file_output_compression, by_stream_gz)
{
    std::ostringstream out;
    {
        seqan3::contrib::gz_ostream compout{out};
        compression_by_stream_impl(compout);
    }
    std::string buffer = out.str();
    buffer[9] = '\x00'; // zero out OS byte
    SEQAN3_TEST_GTEST_SKIP_ZLIB_DEFLATE;
    EXPECT_EQ(buffer, expected_gz);
}

TEST_F(structure_file_output_compression, by_filename_bgzf)
{
    seqan3::test::tmp_directory tmp;
    auto filename = tmp.path() / "structure_file_output_test.dbn.bgzf";
    std::string buffer = compression_by_filename_impl(filename);
    buffer[9] = '\x00'; // zero out OS byte
    SEQAN3_TEST_GTEST_SKIP_ZLIB_DEFLATE;
    EXPECT_EQ(buffer, expected_bgzf);
}

TEST_F(structure_file_output_compression, by_stream_bgzf)
{
    std::ostringstream out;
    {
        seqan3::contrib::bgzf_ostream compout{out};
        compression_by_stream_impl(compout);
    }
    std::string buffer = out.str();
    buffer[9] = '\x00'; // zero out OS byte
    SEQAN3_TEST_GTEST_SKIP_ZLIB_DEFLATE;
    EXPECT_EQ(buffer, expected_bgzf);
}
#endif

#if defined(SEQAN3_HAS_BZIP2)
std::string expected_bz2{
    '\x42', '\x5A', '\x68', '\x39', '\x31', '\x41', '\x59', '\x26', '\x53', '\x59', '\xC7', '\x0B', '\xB5', '\x7F',
    '\x00', '\x00', '\x36', '\x5F', '\x80', '\x6E', '\x50', '\x40', '\x63', '\xEC', '\x81', '\x2A', '\xC3', '\x5A',
    '\x00', '\xAA', '\x26', '\x5D', '\x40', '\x30', '\x00', '\xB8', '\x84', '\x53', '\xC5', '\x00', '\x68', '\x00',
    '\x03', '\x40', '\x34', '\x69', '\xEA', '\x18', '\x01', '\x93', '\x4D', '\x06', '\x43', '\x04', '\x34', '\xC4',
    '\x68', '\xC0', '\x94', '\xD3', '\x52', '\xA7', '\xEA', '\x9B', '\x14', '\xF2', '\x69', '\x1E', '\xA7', '\xA9',
    '\xEA', '\x68', '\xC0', '\x23', '\x4D', '\x35', '\x85', '\x85', '\xCA', '\x54', '\xA4', '\x4F', '\xB6', '\x4C',
    '\xD9', '\xCB', '\x3C', '\xCD', '\x51', '\x11', '\xE5', '\x16', '\xEB', '\x96', '\x5A', '\x11', '\x7E', '\x14',
    '\xC1', '\x50', '\xCB', '\x07', '\x06', '\x2B', '\x15', '\x01', '\x5B', '\x6E', '\xD5', '\x48', '\x26', '\xEA',
    '\xCA', '\x37', '\x7B', '\xE7', '\xE9', '\x9E', '\xDD', '\x0D', '\x2B', '\x79', '\xF1', '\xF4', '\xB6', '\x8B',
    '\x78', '\xB2', '\x4D', '\x0A', '\x53', '\x43', '\x4D', '\x0D', '\x48', '\xD0', '\x98', '\xDC', '\xC4', '\xC4',
    '\x8C', '\x7F', '\x69', '\x94', '\x48', '\xA2', '\x99', '\x15', '\x53', '\xA1', '\x44', '\xC1', '\x31', '\x02',
    '\x5A', '\xF5', '\x91', '\xA7', '\x00', '\x40', '\x10', '\xC2', '\x66', '\x06', '\x02', '\xE0', '\x81', '\x10',
    '\x09', '\x94', '\x46', '\x6E', '\x8E', '\xBD', '\x26', '\x2C', '\xED', '\x8D', '\x97', '\xE4', '\x47', '\xD1',
    '\x4A', '\x42', '\x0F', '\xC5', '\xDC', '\x91', '\x4E', '\x14', '\x24', '\x31', '\xC2', '\xED', '\x5F', '\xC0'};

TEST_F(structure_file_output_compression, by_filename_bz2)
{
    seqan3::test::tmp_directory tmp;
    auto filename = tmp.path() / "structure_file_output_test.dbn.bz2";
    std::string buffer = compression_by_filename_impl(filename);
    EXPECT_EQ(buffer, expected_bz2);
}

TEST_F(structure_file_output_compression, by_stream_bz2)
{
    std::ostringstream out;
    {
        seqan3::contrib::bz2_ostream compout{out};
        compression_by_stream_impl(compout);
    }
    EXPECT_EQ(out.str(), expected_bz2);
}
#endif
