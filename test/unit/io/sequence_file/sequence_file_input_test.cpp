// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/convert.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <seqan3/test/tmp_filename.hpp>

using namespace seqan3;

TEST(sequence_file_input_iterator, concepts)
{
    using it_t = typename sequence_file_input<>::iterator;
    using sen_t = typename sequence_file_input<>::sentinel;

    EXPECT_TRUE((std::input_iterator<it_t>));
    EXPECT_TRUE((std::sentinel_for<sen_t, it_t>));
}

struct sequence_file_input_f : public ::testing::Test
{
    std::string input
    {
        "> TEST 1\n"
        "ACGT\n"
        ">Test2\n"
        "AGGCTGN\n"
        "> Test3\n"
        "GGAGTATAATATATATATATATAT\n"
    };

    dna5_vector seq_comp[3]
    {
        "ACGT"_dna5,
        "AGGCTGN"_dna5,
        "GGAGTATAATATATATATATATAT"_dna5
    };

    std::string id_comp[3]
    {
        "TEST 1",
        "Test2",
        "Test3"
    };
};

TEST_F(sequence_file_input_f, concepts)
{
    using t = sequence_file_input<>;
    EXPECT_TRUE((std::ranges::input_range<t>));

    using ct = sequence_file_input<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::input_range<ct>));
}

TEST_F(sequence_file_input_f, construct_by_filename)
{
    /* just the filename */
    {
        test::tmp_filename filename{"sequence_file_input_constructor.fasta"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        }

        EXPECT_NO_THROW( sequence_file_input<>{filename.get_path()} );
    }

    // correct format check is done by tests of that format

    /* wrong extension */
    {
        test::tmp_filename filename{"sequence_file_input_constructor.xyz"};
        std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        EXPECT_THROW( sequence_file_input<>{filename.get_path()} ,
                      unhandled_extension_error );
    }

    /* non-existent file */
    {
        EXPECT_THROW( sequence_file_input<>{"/dev/nonexistant/foobarOOO"},
                      file_open_error);
    }

    /* filename + fields */
    {
        test::tmp_filename filename{"sequence_file_input_constructor.fasta"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        }

        EXPECT_NO_THROW(( sequence_file_input<sequence_file_input_default_traits_dna,
                                             fields<field::seq>,
                                             type_list<format_fasta>>{filename.get_path(), fields<field::seq>{}} ));
    }
}

TEST_F(sequence_file_input_f, construct_from_stream)
{
    /* stream + format_tag */
    EXPECT_NO_THROW(( sequence_file_input<sequence_file_input_default_traits_dna,
                                         fields<field::seq, field::id, field::qual>,
                                         type_list<format_fasta>>{std::istringstream{input},
                                                                  format_fasta{}} ));

    /* stream + format_tag + fields */
    EXPECT_NO_THROW(( sequence_file_input<sequence_file_input_default_traits_dna,
                                          fields<field::seq, field::id, field::qual>,
                                          type_list<format_fasta>>{std::istringstream{input},
                                                                   format_fasta{},
                                                                   fields<field::seq, field::id, field::qual>{}} ));
}

TEST_F(sequence_file_input_f, default_template_args_and_deduction_guides)
{
    using comp0 = sequence_file_input_default_traits_dna;
    using comp1 = fields<field::seq, field::id, field::qual>;
    using comp2 = type_list<format_embl, format_fasta, format_fastq, format_genbank, format_sam>;
    using comp3 = char;

    /* default template args */
    {
        using t = sequence_file_input<>;
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided filename constructor */
    {
        test::tmp_filename filename{"sequence_file_input_constructor.fasta"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        }

        sequence_file_input fin{filename.get_path()};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided filename constructor + custom fields */
    {
        test::tmp_filename filename{"sequence_file_input_constructor.fasta"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        }

        sequence_file_input fin{filename.get_path(), fields<field::seq>{}};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::seq>>));                   // changed
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided stream constructor */
    {
        std::istringstream ext{input};
        sequence_file_input fin{ext, format_fasta{}};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<format_fasta>>));// changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided stream temporary constructor */
    {
        sequence_file_input fin{std::istringstream{input}, format_fasta{}};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<format_fasta>>));// changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }
}

TEST_F(sequence_file_input_f, empty_file)
{
    test::tmp_filename filename{"empty.fasta"};
    std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};

    sequence_file_input fin{filename.get_path()};

    EXPECT_EQ(fin.begin(), fin.end());
}

TEST_F(sequence_file_input_f, empty_stream)
{
    sequence_file_input fin{std::istringstream{std::string{}}, format_fasta{}};

    EXPECT_EQ(fin.begin(), fin.end());
}

TEST_F(sequence_file_input_f, record_reading)
{
    /* record based reading */
    sequence_file_input fin{std::istringstream{input}, format_fasta{}};

    size_t counter = 0;
    for (auto & rec : fin)
    {
        EXPECT_TRUE((std::ranges::equal(get<field::seq>(rec), seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(get<field::id>(rec),  id_comp[counter])));
        EXPECT_TRUE(empty(get<field::qual>(rec)));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

TEST_F(sequence_file_input_f, record_reading_struct_bind)
{
    /* record based reading */
    sequence_file_input fin{std::istringstream{input}, format_fasta{}};

    size_t counter = 0;
    for (auto & [ seq, id, qual ] : fin)
    {
        EXPECT_TRUE((std::ranges::equal(seq, seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(id,  id_comp[counter])));
        EXPECT_TRUE(empty(qual));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

TEST_F(sequence_file_input_f, record_reading_custom_fields)
{
    /* record based reading */
    sequence_file_input fin{std::istringstream{input},
                         format_fasta{},
                         fields<field::id, field::seq_qual>{}};

    size_t counter = 0;
    for (auto & [ id, seq_qual ] : fin)
    {
        EXPECT_TRUE((std::ranges::equal(seq_qual | views::convert<dna5>, seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(id,  id_comp[counter])));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

TEST_F(sequence_file_input_f, record_reading_custom_options)
{
    std::istringstream istream{std::string
    {
        "> ID1 lala\n"
        "ACGTTTTTTTTTTTTTTT\n"
        "> ID2\n"
        "ACGTTTTTTT\n"
        "> ID3 lala\n"
        "ACGTTTA\n"
    }};

    /* record based reading */
    sequence_file_input fin{istream, format_fasta{}};
    fin.options.truncate_ids = true;

    auto it = fin.begin();
    EXPECT_EQ(get<field::id>(*it), "ID1");
    ++it;
    EXPECT_EQ(get<field::id>(*it), "ID2");
    ++it;
    EXPECT_EQ(get<field::id>(*it), "ID3");
}

TEST_F(sequence_file_input_f, file_view)
{
    sequence_file_input fin{std::istringstream{input}, format_fasta{}};

    auto minimum_length_filter = std::views::filter([] (auto const & rec)
    {
        return size(get<field::seq>(rec)) >= 5;
    });

    size_t counter = 1; // the first record will be filtered out
    for (auto & rec : fin | minimum_length_filter)
    {
        EXPECT_TRUE((std::ranges::equal(get<field::seq>(rec), seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(get<field::id>(rec),  id_comp[counter])));
        EXPECT_TRUE(empty(get<field::qual>(rec)));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

// ----------------------------------------------------------------------------
// decompression
// ----------------------------------------------------------------------------

template <typename fixture_t, typename input_file_t>
void decompression_impl(fixture_t & fix, input_file_t & fin)
{
    size_t counter = 0;
    for (auto & rec : fin)
    {
        EXPECT_TRUE((std::ranges::equal(get<field::seq>(rec), fix.seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(get<field::id>(rec),  fix.id_comp[counter])));
        EXPECT_TRUE(empty(get<field::qual>(rec)));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

#ifdef SEQAN3_HAS_ZLIB
std::string input_gz
{
    '\x1F','\x8B','\x08','\x00','\x33','\xBF','\x13','\x5C','\x00','\x03','\xB3','\x53','\x08','\x71','\x0D','\x0E',
    '\x51','\x30','\xE4','\x72','\x74','\x76','\x0F','\xE1','\xB2','\x0B','\x49','\x2D','\x2E','\x31','\xE2','\x72',
    '\x74','\x77','\x77','\x0E','\x71','\xF7','\xE3','\xB2','\x53','\x00','\xF1','\x8D','\xB9','\xDC','\xDD','\x1D',
    '\xDD','\x43','\x1C','\x43','\x1C','\x1D','\x43','\x50','\x21','\x17','\x00','\xEF','\x24','\xC2','\xE9','\x3E',
    '\x00','\x00','\x00',
};

TEST_F(sequence_file_input_f, decompression_by_filename_gz)
{
    test::tmp_filename filename{"sequence_file_output_test.fasta.gz"};

    {
        std::ofstream of{filename.get_path(), std::ios::binary};

        std::copy(begin(input_gz), end(input_gz), std::ostreambuf_iterator<char>{of});
    }

    sequence_file_input fin{filename.get_path()};

    decompression_impl(*this, fin);
}

TEST_F(sequence_file_input_f, decompression_by_stream_gz)
{
    sequence_file_input fin{std::istringstream{input_gz}, format_fasta{}};

    decompression_impl(*this, fin);
}

TEST_F(sequence_file_input_f, read_empty_gz_file)
{
    std::string empty_zipped_file
    {
        '\x1f', '\x8b', '\x08', '\x08', '\x5a', '\x07', '\x98', '\x5c',
        '\x00', '\x03', '\x66', '\x6f', '\x6f', '\x00', '\x03', '\x00',
        '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00'
    };
    sequence_file_input fin{std::istringstream{empty_zipped_file}, format_fasta{}};

    EXPECT_TRUE(fin.begin() == fin.end());
}

std::string input_bgzf
{
    '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43',
    '\x02', '\x00', '\x4A', '\x00', '\xB3', '\x53', '\x08', '\x71', '\x0D', '\x0E', '\x51', '\x30', '\xE4', '\x72',
    '\x74', '\x76', '\x0F', '\xE1', '\xB2', '\x0B', '\x49', '\x2D', '\x2E', '\x31', '\xE2', '\x72', '\x74', '\x77',
    '\x77', '\x0E', '\x71', '\xF7', '\xE3', '\xB2', '\x53', '\x00', '\xF1', '\x8D', '\xB9', '\xDC', '\xDD', '\x1D',
    '\xDD', '\x43', '\x1C', '\x43', '\x1C', '\x1D', '\x43', '\x50', '\x21', '\x17', '\x00', '\xEF', '\x24', '\xC2',
    '\xE9', '\x3E', '\x00', '\x00', '\x00', '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00',
    '\xFF', '\x06', '\x00', '\x42', '\x43', '\x02', '\x00', '\x1B', '\x00', '\x03', '\x00', '\x00', '\x00', '\x00',
    '\x00', '\x00', '\x00', '\x00', '\x00'
};

TEST_F(sequence_file_input_f, bgzf_decompression_by_filename_bgzf)
{
    test::tmp_filename filename{"sequence_file_output_test.fasta.bgzf"};

    {
        std::ofstream of{filename.get_path(), std::ios::binary};

        std::copy(input_bgzf.begin(), input_bgzf.end(), std::ostreambuf_iterator<char>{of});
    }

    sequence_file_input fin{filename.get_path()};

    decompression_impl(*this, fin);
}

TEST_F(sequence_file_input_f, bgzf_decompression_by_filename_gz)
{
    test::tmp_filename filename{"sequence_file_output_test.fasta.gz"};

    {
        std::ofstream of{filename.get_path(), std::ios::binary};
        std::copy(input_bgzf.begin(), input_bgzf.end(), std::ostreambuf_iterator<char>{of});
    }

    sequence_file_input fin{filename.get_path()};
    decompression_impl(*this, fin);
}

TEST_F(sequence_file_input_f, decompression_by_stream_bgzf)
{
    sequence_file_input fin{std::istringstream{input_bgzf}, format_fasta{}};

    decompression_impl(*this, fin);
}

TEST_F(sequence_file_input_f, read_empty_bgzf_file)
{
    std::string empty_bgzf_file
    {
        '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF',
        '\x06', '\x00', '\x42', '\x43', '\x02', '\x00', '\x1B', '\x00', '\x03', '\x00',
        '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00',
    };
    sequence_file_input fin{std::istringstream{empty_bgzf_file}, format_fasta{}};

    EXPECT_TRUE(fin.begin() == fin.end());
}
#endif

#ifdef SEQAN3_HAS_BZIP2
std::string input_bz2
{
    '\x42','\x5A','\x68','\x39','\x31','\x41','\x59','\x26','\x53','\x59','\x8D','\xD7','\xE7','\xD6','\x00','\x00',
    '\x06','\x5F','\x80','\x00','\x10','\x40','\x00','\x38','\x01','\x2A','\x81','\x0C','\x00','\x02','\x00','\x0C',
    '\x00','\x20','\x00','\x54','\x44','\x34','\xC0','\x00','\x4A','\x9B','\x44','\x68','\x9E','\x48','\x5D','\x34',
    '\x67','\x4F','\x24','\xFC','\x6F','\x10','\xC5','\xA0','\x3C','\x12','\x61','\xDD','\xE9','\x45','\xA5','\xD4',
    '\x26','\x31','\xBC','\xF1','\x49','\x61','\x81','\xA2','\xEE','\x48','\xA7','\x0A','\x12','\x11','\xBA','\xFC',
    '\xFA','\xC0'
};

TEST_F(sequence_file_input_f, decompression_by_filename_bz2)
{
    test::tmp_filename filename{"sequence_file_output_test.fasta.bz2"};

    {
        std::ofstream of{filename.get_path(), std::ios::binary};

        std::copy(begin(input_bz2), end(input_bz2), std::ostreambuf_iterator<char>{of});
    }

    sequence_file_input fin{filename.get_path()};

    decompression_impl(*this, fin);
}

TEST_F(sequence_file_input_f, decompression_by_stream_bz2)
{
    sequence_file_input fin{std::istringstream{input_bz2}, format_fasta{}};

    decompression_impl(*this, fin);
}

TEST_F(sequence_file_input_f, read_empty_bz2_file)
{
    std::string empty_zipped_file
    {
        '\x42', '\x5a', '\x68', '\x39', '\x17', '\x72', '\x45', '\x38', '\x50', '\x90', '\x00', '\x00', '\x00', '\x00'
    };
    sequence_file_input fin{std::istringstream{empty_zipped_file}, format_fasta{}};

    EXPECT_TRUE(fin.begin() == fin.end());
}
#endif
