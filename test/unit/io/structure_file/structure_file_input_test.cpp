// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <fstream>
#include <iterator>
#include <ranges>
#include <sstream>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/tmp_directory.hpp>
#include <seqan3/utility/views/convert.hpp>

using seqan3::operator""_rna5;
using seqan3::operator""_wuss51;

TEST(general, concepts)
{
    using it_t = typename seqan3::structure_file_input<>::iterator;
    using sen_t = typename seqan3::structure_file_input<>::sentinel;

    EXPECT_TRUE((std::input_iterator<it_t>));
    EXPECT_TRUE((std::sentinel_for<sen_t, it_t>));
}

struct structure_file_input_class : public ::testing::Test
{
    using comp0 = seqan3::structure_file_input_default_traits_rna;
    using comp1 = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>;
    using comp2 = seqan3::type_list<seqan3::format_vienna>;
    using comp3 = char;

    seqan3::test::tmp_directory directory_tmp{};

    seqan3::test::sandboxed_path create_file(char const * contents)
    {
        auto filename = directory_tmp.path() / "structure_file_input_constructor.dbn";
        std::ofstream filecreator{filename, std::ios::out | std::ios::binary};
        filecreator << contents; // must contain at least one record
        return filename;
    }
};

TEST_F(structure_file_input_class, concepts)
{
    using t = seqan3::structure_file_input<>;
    EXPECT_TRUE((std::ranges::input_range<t>));

    using ct = seqan3::structure_file_input<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::input_range<ct>));
}

TEST_F(structure_file_input_class, construct_by_filename)
{
    /* just the filename */
    {
        seqan3::test::tmp_directory tmp{};
        auto filename = tmp.path() / "structure_file_input_constructor.dbn";

        {
            std::ofstream filecreator{filename, std::ios::out | std::ios::binary};
        }

        EXPECT_NO_THROW(seqan3::structure_file_input<>{filename});
    }

    // correct format check is done by tests of that format

    /* wrong extension */
    {
        seqan3::test::tmp_directory tmp{};
        auto filename = tmp.path() / "structure_file_input_constructor.xyz";

        std::ofstream filecreator{filename, std::ios::out | std::ios::binary};
        EXPECT_THROW(seqan3::structure_file_input<>{filename}, seqan3::unhandled_extension_error);
    }

    /* non-existent file*/
    {
        EXPECT_THROW(seqan3::structure_file_input<>{"/dev/nonexistant/foobarOOO"}, seqan3::file_open_error);
    }

    /* filename + fields */
    {
        using fields_seq = seqan3::fields<seqan3::field::seq>;

        auto filename = create_file("> ID\nACGU\n....\n");
        EXPECT_NO_THROW(
            (seqan3::structure_file_input<seqan3::structure_file_input_default_traits_rna,
                                          fields_seq,
                                          seqan3::type_list<seqan3::format_vienna>>{filename, fields_seq{}}));
    }
}

TEST_F(structure_file_input_class, construct_from_stream)
{
    using fields_seq_id_structure = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::structure>;
    std::string str = "> ID\nACGU\n....\n";

    /* stream + format_tag */
    EXPECT_NO_THROW((seqan3::structure_file_input<seqan3::structure_file_input_default_traits_rna,
                                                  fields_seq_id_structure,
                                                  seqan3::type_list<seqan3::format_vienna>>{std::istringstream{str},
                                                                                            seqan3::format_vienna{}}));

    /* stream + format_tag + fields */
    EXPECT_NO_THROW(
        (seqan3::structure_file_input<seqan3::structure_file_input_default_traits_rna,
                                      fields_seq_id_structure,
                                      seqan3::type_list<seqan3::format_vienna>>{std::istringstream{str},
                                                                                seqan3::format_vienna{},
                                                                                fields_seq_id_structure{}}));
}

TEST_F(structure_file_input_class, default_template_args)
{
    /* default template args */
    using t = seqan3::structure_file_input<>;
    EXPECT_TRUE((std::is_same_v<typename t::traits_type, comp0>));
    EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
    EXPECT_TRUE((std::is_same_v<typename t::valid_formats, comp2>));
    EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
}

TEST_F(structure_file_input_class, guided_filename_constructor)
{
    /* guided filename constructor */
    auto filename = create_file("> ID\nACGU\n....\n");
    seqan3::structure_file_input fin{filename};

    using t = decltype(fin);
    EXPECT_TRUE((std::is_same_v<typename t::traits_type, comp0>));
    EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
    EXPECT_TRUE((std::is_same_v<typename t::valid_formats, comp2>));
    EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
}

TEST_F(structure_file_input_class, guided_filename_constructor_and_custom_fields)
{
    /* guided filename constructor + custom fields */
    auto filename = create_file("> ID\nACGU\n....\n");
    seqan3::structure_file_input fin{filename, seqan3::fields<seqan3::field::seq>{}};

    using t = decltype(fin);
    EXPECT_TRUE((std::is_same_v<typename t::traits_type, comp0>));
    EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, seqan3::fields<seqan3::field::seq>>));
    EXPECT_TRUE((std::is_same_v<typename t::valid_formats, comp2>));
    EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
}

TEST_F(structure_file_input_class, guided_stream_constructor)
{
    /* guided stream constructor */
    seqan3::structure_file_input fin{std::istringstream{"> ID\nACGU\n....\n"}, seqan3::format_vienna{}};

    using t = decltype(fin);
    EXPECT_TRUE((std::is_same_v<typename t::traits_type, comp0>));
    EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
    EXPECT_TRUE((std::is_same_v<typename t::valid_formats, seqan3::type_list<seqan3::format_vienna>>));
    EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
}

TEST_F(structure_file_input_class, amino_acids_traits)
{
    auto filename = create_file("> ID\nACEW\nHHHH\n");
    seqan3::structure_file_input<seqan3::structure_file_input_default_traits_aa> fin{filename};

    using t = decltype(fin);
    EXPECT_TRUE((std::is_same_v<typename t::traits_type, seqan3::structure_file_input_default_traits_aa>));
    EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
    EXPECT_TRUE((std::is_same_v<typename t::valid_formats, comp2>));
    EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
}

TEST_F(structure_file_input_class, modified_traits)
{
    auto filename = create_file("> ID\nACGU\n....\n");
    //! [structure_file_input_class mod_traits]
    struct my_traits : seqan3::structure_file_input_default_traits_rna
    {
        using seq_alphabet = seqan3::rna4; // instead of rna5
    };

    seqan3::structure_file_input<my_traits> fin{filename};
    //! [structure_file_input_class mod_traits]

    using t = decltype(fin);
    EXPECT_TRUE((std::is_same_v<typename t::traits_type, my_traits>));
    EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
    EXPECT_TRUE((std::is_same_v<typename t::valid_formats, comp2>));
    EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
}

struct structure_file_input_read : public ::testing::Test
{
    size_t const num_records = 2ul;

    std::string const input{">S.cerevisiae_tRNA-PHE M10740/1-73\n"
                            "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
                            "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
                            "> example\n"
                            "UUGGAGUACACAACCUGUACACUCUUUC\n"
                            "..(((((..(((...)))..)))))... (-3.71)\n"};

    seqan3::rna5_vector const seq_comp[2]{
        "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA"_rna5,
        "UUGGAGUACACAACCUGUACACUCUUUC"_rna5};

    std::string const id_comp[2]{"S.cerevisiae_tRNA-PHE M10740/1-73", "example"};

    double const energy_comp[2]{-17.5, -3.71};

    std::vector<seqan3::wuss51> const structure_comp[2]{
        "(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))."_wuss51,
        "..(((((..(((...)))..)))))..."_wuss51};

    std::vector<uint8_t> const interaction_comp[2]{{71, 70, 69, 68, 67, 66, 65, 24, 23, 22, 21, 12, 11, 10,
                                                    9,  42, 41, 40, 39, 29, 28, 27, 26, 64, 63, 62, 61, 60,
                                                    52, 51, 50, 49, 48, 6,  5,  4,  3,  2,  1,  0},
                                                   {24, 23, 22, 21, 20, 17, 16, 15, 11, 10, 9, 6, 5, 4, 3, 2}};

    template <typename bpp_type>
    void bpp_test(bpp_type & bpp, std::vector<uint8_t> const & bpp_comp)
    {
        size_t idx = 0ul;
        auto interactions = bpp
                          | std::views::filter(
                                [](auto & set)
                                {
                                    return set.size() == 1;
                                });
        for (auto & elem : interactions)
        {
            EXPECT_EQ(elem.begin()->second, bpp_comp[idx++]);
        }
        EXPECT_EQ(idx, bpp_comp.size());
    }

#if defined(SEQAN3_HAS_ZLIB) || defined(SEQAN3_HAS_BZIP2)
    template <typename input_file_t>
    void decompression_impl(input_file_t & fin)
    {
        size_t counter = 0ul;
        for (auto & rec : fin)
        {
            EXPECT_RANGE_EQ(rec.sequence(), seq_comp[counter]);
            EXPECT_RANGE_EQ(rec.id(), id_comp[counter]);
            EXPECT_RANGE_EQ(rec.sequence_structure(), structure_comp[counter]);

            ++counter;
        }

        EXPECT_EQ(counter, num_records);
    }
#endif
};

TEST_F(structure_file_input_read, empty_file)
{
    seqan3::test::tmp_directory tmp{};
    auto filename = tmp.path() / "empty.dbn";
    std::ofstream filecreator{filename, std::ios::out | std::ios::binary};

    seqan3::structure_file_input fin{filename};

    EXPECT_EQ(fin.begin(), fin.end());
}

TEST_F(structure_file_input_read, empty_stream)
{
    seqan3::structure_file_input fin{std::istringstream{std::string{}}, seqan3::format_vienna{}};

    EXPECT_EQ(fin.begin(), fin.end());
}

TEST_F(structure_file_input_read, record_general)
{
    /* record based reading */
    seqan3::structure_file_input fin{std::istringstream{input}, seqan3::format_vienna{}};

    size_t counter = 0ul;
    for (auto & rec : fin)
    {
        EXPECT_RANGE_EQ(rec.sequence(), seq_comp[counter]);
        EXPECT_RANGE_EQ(rec.id(), id_comp[counter]);
        EXPECT_RANGE_EQ(rec.sequence_structure(), structure_comp[counter]);

        ++counter;
    }
    EXPECT_EQ(counter, num_records);
}

TEST_F(structure_file_input_read, record_struct_bind)
{
    /* record based reading */
    seqan3::structure_file_input fin{std::istringstream{input},
                                     seqan3::format_vienna{},
                                     seqan3::fields<seqan3::field::seq,
                                                    seqan3::field::id,
                                                    seqan3::field::bpp,
                                                    seqan3::field::structure,
                                                    seqan3::field::energy>{}};

    size_t counter = 0ul;
    for (auto & [sequence, id, bpp, structure, energy] : fin)
    {
        EXPECT_RANGE_EQ(sequence, seq_comp[counter]);
        EXPECT_RANGE_EQ(id, id_comp[counter]);
        EXPECT_RANGE_EQ(structure, structure_comp[counter]);
        EXPECT_DOUBLE_EQ(energy.value(), energy_comp[counter]);
        bpp_test(bpp, interaction_comp[counter]);
        ++counter;
    }
    EXPECT_EQ(counter, num_records);
}

TEST_F(structure_file_input_read, record_custom_fields)
{
    /* record based reading */
    seqan3::structure_file_input fin{std::istringstream{input},
                                     seqan3::format_vienna{},
                                     seqan3::fields<seqan3::field::id, seqan3::field::structured_seq>{}};

    size_t counter = 0ul;
    for (auto & [id, seq_structure] : fin)
    {
        EXPECT_RANGE_EQ(id, id_comp[counter]);
        EXPECT_RANGE_EQ(seq_structure | seqan3::views::convert<seqan3::rna5>, seq_comp[counter]);
        EXPECT_RANGE_EQ(seq_structure | seqan3::views::convert<seqan3::wuss51>, structure_comp[counter]);
        ++counter;
    }
    EXPECT_EQ(counter, num_records);
}

TEST_F(structure_file_input_read, record_file_view)
{
    seqan3::structure_file_input fin{std::istringstream{input},
                                     seqan3::format_vienna{},
                                     seqan3::fields<seqan3::field::seq,
                                                    seqan3::field::id,
                                                    seqan3::field::bpp,
                                                    seqan3::field::structure,
                                                    seqan3::field::energy>{}};

    auto minimum_length_filter = std::views::filter(
        [](auto const & rec)
        {
            return size(rec.sequence()) >= 5;
        });

    size_t counter = 0ul; // the first record will be filtered out
    for (auto & rec : fin | minimum_length_filter)
    {
        EXPECT_RANGE_EQ(rec.sequence(), seq_comp[counter]);
        EXPECT_RANGE_EQ(rec.id(), id_comp[counter]);
        bpp_test(rec.base_pair_probability_matrix(), interaction_comp[counter]);
        EXPECT_RANGE_EQ(rec.sequence_structure(), structure_comp[counter]);
        EXPECT_DOUBLE_EQ(rec.energy().value(), energy_comp[counter]);

        ++counter;
    }
    EXPECT_EQ(counter, num_records);
}

// ----------------------------------------------------------------------------
// decompression
// ----------------------------------------------------------------------------

#if defined(SEQAN3_HAS_ZLIB)
std::string input_gz{
    '\x1F', '\x8B', '\x08', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x03', '\x55', '\x8E', '\xC1', '\x0A',
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

TEST_F(structure_file_input_read, decompression_by_filename_gz)
{
    seqan3::test::tmp_directory tmp{};
    auto filename = tmp.path() / "structure_file_input_test.dbn.gz";

    {
        std::ofstream of{filename, std::ios::binary};

        std::copy(input_gz.begin(), input_gz.end(), std::ostreambuf_iterator<char>{of});
    }

    seqan3::structure_file_input fin{filename};

    decompression_impl(fin);
}

TEST_F(structure_file_input_read, decompression_by_stream_gz)
{
    seqan3::structure_file_input fin{std::istringstream{input_gz}, seqan3::format_vienna{}};

    decompression_impl(fin);
}

TEST_F(structure_file_input_read, read_empty_gz_file)
{
    std::string empty_zipped_file{'\x1f', '\x8b', '\x08', '\x08', '\x5a', '\x07', '\x98', '\x5c',
                                  '\x00', '\x03', '\x66', '\x6f', '\x6f', '\x00', '\x03', '\x00',
                                  '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00'};
    seqan3::structure_file_input fin{std::istringstream{empty_zipped_file}, seqan3::format_vienna{}};

    EXPECT_TRUE(fin.begin() == fin.end());
}

std::string input_bgzf{
    '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43',
    '\x02', '\x00', '\xB6', '\x00', '\x55', '\x8E', '\xC1', '\x0A', '\x02', '\x31', '\x0C', '\x44', '\xEF', '\xF9',
    '\x8A', '\x1E', '\xDB', '\xC3', '\xD6', '\x96', '\x55', '\x7A', '\x5B', '\x08', '\x41', '\xE2', '\x45', '\x11',
    '\x25', '\x67', '\x59', '\xA4', '\x87', '\x05', '\x05', '\x51', '\x11', '\x3F', '\xDF', '\xB4', '\x15', '\x59',
    '\x27', '\x50', '\x3A', '\x8F', '\x30', '\x93', '\xE1', '\xE8', '\xCF', '\xF9', '\x9E', '\x5F', '\xD3', '\x63',
    '\x1A', '\xF3', '\xE9', '\x79', '\xD8', '\x61', '\xB7', '\xDF', '\xAC', '\xCD', '\x36', '\x86', '\xB4', '\x0C',
    '\x8B', '\xD8', '\xA5', '\x1E', '\x98', '\x98', '\x51', '\x44', '\x90', '\x49', '\x08', '\x59', '\x84', '\xD5',
    '\xEB', '\x10', '\x93', '\x5A', '\x24', '\x61', '\xC4', '\xBA', '\x50', '\x30', '\x0B', '\x29', '\xD0', '\x11',
    '\x52', '\xA6', '\x0B', '\x65', '\x05', '\x8B', '\x23', '\x04', '\xDB', '\xE4', '\x7D', '\x7B', '\x9B', '\x9C',
    '\xEA', '\x0F', '\x34', '\x52', '\x64', '\x67', '\xD8', '\xCD', '\xE4', '\x8D', '\xED', '\x62', '\xF2', '\xAB',
    '\xE0', '\x60', '\x30', '\xF9', '\x3D', '\x5E', '\x6F', '\x97', '\x0C', '\xAD', '\x5F', '\xB0', '\x54', '\x62',
    '\x3D', '\xA2', '\x7C', '\x85', '\xF4', '\x30', '\x82', '\x5F', '\x56', '\x0B', '\xAC', '\x05', '\xEE', '\xDB',
    '\xA3', '\x61', '\xBD', '\x4F', '\xD1', '\xC1', '\x07', '\x38', '\xAB', '\x49', '\x82', '\x0C', '\x01', '\x00',
    '\x00', '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42',
    '\x43', '\x02', '\x00', '\x1B', '\x00', '\x03', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00',
    '\x00'};

TEST_F(structure_file_input_read, decompression_by_filename_bgzf)
{
    seqan3::test::tmp_directory tmp{};
    auto filename = tmp.path() / "structure_file_input_test.dbn.bgzf";

    {
        std::ofstream of{filename, std::ios::binary};

        std::copy(input_gz.begin(), input_gz.end(), std::ostreambuf_iterator<char>{of});
    }

    seqan3::structure_file_input fin{filename};

    decompression_impl(fin);
}

TEST_F(structure_file_input_read, decompression_by_stream_bgzf)
{
    seqan3::structure_file_input fin{std::istringstream{input_gz}, seqan3::format_vienna{}};

    decompression_impl(fin);
}

TEST_F(structure_file_input_read, read_empty_bgzf_file)
{
    std::string empty_bgzf_file{'\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF',
                                '\x06', '\x00', '\x42', '\x43', '\x02', '\x00', '\x1B', '\x00', '\x03', '\x00',
                                '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00'};
    seqan3::structure_file_input fin{std::istringstream{empty_bgzf_file}, seqan3::format_vienna{}};

    EXPECT_TRUE(fin.begin() == fin.end());
}

#endif

#if defined(SEQAN3_HAS_BZIP2)
std::string input_bz2{
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

TEST_F(structure_file_input_read, decompression_by_filename_bz2)
{
    seqan3::test::tmp_directory tmp{};
    auto filename = tmp.path() / "structure_file_input_test.dbn.bz2";

    {
        std::ofstream of{filename, std::ios::binary};

        std::copy(input_bz2.begin(), input_bz2.end(), std::ostreambuf_iterator<char>{of});
    }

    seqan3::structure_file_input fin{filename};

    decompression_impl(fin);
}

TEST_F(structure_file_input_read, decompression_by_stream_bz2)
{
    seqan3::structure_file_input fin{std::istringstream{input_bz2}, seqan3::format_vienna{}};

    decompression_impl(fin);
}

TEST_F(structure_file_input_read, read_empty_bz2_file)
{
    std::string empty_zipped_file{'\x42',
                                  '\x5a',
                                  '\x68',
                                  '\x39',
                                  '\x17',
                                  '\x72',
                                  '\x45',
                                  '\x38',
                                  '\x50',
                                  '\x90',
                                  '\x00',
                                  '\x00',
                                  '\x00',
                                  '\x00'};
    seqan3::structure_file_input fin{std::istringstream{empty_zipped_file}, seqan3::format_vienna{}};

    EXPECT_TRUE(fin.begin() == fin.end());
}
#endif
