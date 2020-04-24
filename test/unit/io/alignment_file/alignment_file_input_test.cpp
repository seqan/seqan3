// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/alignment_file/input.hpp>
#include <seqan3/range/views/convert.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <seqan3/test/tmp_filename.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;
using seqan3::operator""_phred42;

using default_fields = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::qual>;

TEST(alignment_file_input_iterator, concepts)
{
    using it_t = typename seqan3::alignment_file_input<>::iterator;
    using sen_t = typename seqan3::alignment_file_input<>::sentinel;

    EXPECT_TRUE((std::input_iterator<it_t>));
    EXPECT_TRUE((std::sentinel_for<sen_t, it_t>));
}

struct alignment_file_input_f : public ::testing::Test
{
    std::string input =
R"(@HD	VN:1.6	SO:unknown	GO:none
@SQ	SN:ref	LN:34
@PG	ID:prog1	PN:cool_program
@CO	This is a comment.
read1	41	ref	1	61	1S1M1D2M	ref	10	300	ACGT	!##$	AS:i:2	NM:i:7
read2	42	ref	2	62	7M1D1M1S	ref	10	300	AGGCTGNAG	!##$&'()*	xy:B:S,3,4,5
read3	43	ref	3	63	1S1M1D4M1D1M1S	ref	10	300	GGAGTATA	!!*+,-./
)";

    std::vector<seqan3::dna5_vector> seq_comp
    {
        "ACGT"_dna5,
        "AGGCTGNAG"_dna5,
        "GGAGTATA"_dna5
    };

    std::vector<std::string> id_comp
    {
        "read1",
        "read2",
        "read3"
    };

    std::vector<std::vector<seqan3::phred42>> qual_comp
    {
        { "!##$"_phred42 },
        { "!##$&'()*"_phred42 },
        { "!!*+,-./"_phred42 },
    };
};

TEST_F(alignment_file_input_f, concepts)
{
    using t = seqan3::alignment_file_input<>;
    EXPECT_TRUE((std::ranges::input_range<t>));

    using ct = seqan3::alignment_file_input<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::input_range<ct>));
}

TEST_F(alignment_file_input_f, construct_by_filename)
{
    /* just the filename */
    {
        seqan3::test::tmp_filename filename{"alignment_file_input_constructor.sam"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        }

        EXPECT_NO_THROW(seqan3::alignment_file_input<>{filename.get_path()} );
    }

    // correct format check is done by tests of that format

    /* wrong extension */
    {
        seqan3::test::tmp_filename filename{"alignment_file_input_constructor.xyz"};
        std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        EXPECT_THROW(seqan3::alignment_file_input<>{filename.get_path()},
                     seqan3::unhandled_extension_error );
    }

    /* non-existent file*/
    {
        EXPECT_THROW(seqan3::alignment_file_input<>{"/dev/nonexistent/foobarOOO"}, seqan3::file_open_error);
    }
    /* non-existent file with reference information*/
    {
        std::vector<std::string> ref_ids{"ref1", "ref2"};
        std::vector<seqan3::dna4_vector> ref_seqs{"ACTG"_dna4, "ACTG"_dna4};
        EXPECT_THROW((seqan3::alignment_file_input{"/dev/nonexistent/foobarOOO", ref_ids, ref_seqs}),
                     seqan3::file_open_error);
    }

    /* filename + fields */
    {
        seqan3::test::tmp_filename filename{"alignment_file_input_constructor.sam"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        }

        using fields_seq = seqan3::fields<seqan3::field::seq>;

        EXPECT_NO_THROW(( seqan3::alignment_file_input<seqan3::alignment_file_input_default_traits<>,
                                                       fields_seq,
                                                       seqan3::type_list<seqan3::format_sam>>{filename.get_path(),
                                                                                              fields_seq{}} ));
    }
}

TEST_F(alignment_file_input_f, construct_from_stream)
{
    /* stream + format_tag */
    EXPECT_NO_THROW(( seqan3::alignment_file_input<seqan3::alignment_file_input_default_traits<>,
                                                   default_fields,
                                                   seqan3::type_list<seqan3::format_sam>>{std::istringstream{input},
                                                                                          seqan3::format_sam{}} ));

    /* stream + format_tag + fields */
    EXPECT_NO_THROW(( seqan3::alignment_file_input<seqan3::alignment_file_input_default_traits<>,
                                                   default_fields,
                                                   seqan3::type_list<seqan3::format_sam>>{std::istringstream{input},
                                                                                          seqan3::format_sam{},
                                                                                          default_fields{}} ));
}

TEST_F(alignment_file_input_f, default_template_args_and_deduction_guides)
{
    using comp0 = seqan3::alignment_file_input_default_traits<>;
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
        using t = seqan3::alignment_file_input<>;
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided filename constructor */
    {
        seqan3::test::tmp_filename filename{"alignment_file_input_constructor.sam"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        }

        seqan3::alignment_file_input fin{filename.get_path()};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided filename constructor + custom fields */
    {
        seqan3::test::tmp_filename filename{"alignment_file_input_constructor.sam"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        }

        seqan3::alignment_file_input fin{filename.get_path(), seqan3::fields<seqan3::field::seq>{}};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, seqan3::fields<seqan3::field::seq>>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided stream constructor */
    {
        std::istringstream ext{input};
        seqan3::alignment_file_input fin{ext, seqan3::format_sam{}};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      seqan3::type_list<seqan3::format_sam>>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }

    /* guided stream temporary constructor */
    {
        seqan3::alignment_file_input fin{std::istringstream{input}, seqan3::format_sam{}};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      seqan3::type_list<seqan3::format_sam>>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type,   comp3>));
    }
}

TEST_F(alignment_file_input_f, empty_file)
{
    seqan3::test::tmp_filename filename{"empty.sam"};
    std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};

    seqan3::alignment_file_input fin{filename.get_path()};

    EXPECT_EQ(fin.begin(), fin.end());
}

TEST_F(alignment_file_input_f, empty_stream)
{
    seqan3::alignment_file_input fin{std::istringstream{std::string{}}, seqan3::format_sam{}};

    EXPECT_EQ(fin.begin(), fin.end());
}

TEST_F(alignment_file_input_f, record_reading)
{
    /* record based reading */
    seqan3::alignment_file_input fin{std::istringstream{input}, seqan3::format_sam{}};

    size_t counter = 0;
    for (auto & rec : fin)
    {
        EXPECT_TRUE((std::ranges::equal(seqan3::get<seqan3::field::seq>(rec), seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(seqan3::get<seqan3::field::id>(rec), id_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(seqan3::get<seqan3::field::qual>(rec), qual_comp[counter])));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

TEST_F(alignment_file_input_f, record_reading_custom_fields)
{
    /* record based reading */
    seqan3::alignment_file_input fin{std::istringstream{input},
                                     seqan3::format_sam{},
                                     seqan3::fields<seqan3::field::id, seqan3::field::seq>{}};

    size_t counter = 0;
    for (auto & [ id, seq ] : fin)
    {
        EXPECT_TRUE((std::ranges::equal(seq, seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(id,  id_comp[counter])));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

TEST_F(alignment_file_input_f, file_view)
{
    seqan3::alignment_file_input fin{std::istringstream{input}, seqan3::format_sam{}};

    auto minimum_length_filter = std::views::filter([] (auto const & rec)
    {
        return size(seqan3::get<seqan3::field::seq>(rec)) >= 5;
    });

    size_t counter = 1; // the first record will be filtered out
    for (auto & rec : fin | minimum_length_filter)
    {
        EXPECT_TRUE((std::ranges::equal(seqan3::get<seqan3::field::seq>(rec), seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(seqan3::get<seqan3::field::id>(rec), id_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(seqan3::get<seqan3::field::qual>(rec), qual_comp[counter])));

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
        EXPECT_TRUE((std::ranges::equal(seqan3::get<seqan3::field::seq>(rec), fix.seq_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(seqan3::get<seqan3::field::id>(rec),  fix.id_comp[counter])));
        EXPECT_TRUE((std::ranges::equal(seqan3::get<seqan3::field::qual>(rec), fix.qual_comp[counter])));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

#ifdef SEQAN3_HAS_ZLIB

std::string input_gz
{
    '\x1F','\x8B','\x08','\x08','\x9D','\x5B','\x38','\x5C','\x00','\x03','\x74','\x65','\x73','\x74','\x2E','\x73',
    '\x61','\x6D','\x00','\x6D','\xCE','\xBF','\x0A','\xC2','\x30','\x10','\xC7','\xF1','\xF9','\xD7','\xB7','\x08',
    '\x15','\xFF','\xD4','\xA8','\xB9','\x24','\xB6','\x90','\x2D','\x56','\xB8','\x29','\x5D','\x92','\x17','\x28',
    '\x58','\xC1','\x35','\x93','\xBE','\xBD','\x21','\xBA','\x08','\x4E','\x07','\xC7','\x7D','\xBE','\x5C','\x5E',
    '\xE6','\x1B','\xC1','\x12','\xF2','\x72','\x07','\xA1','\x27','\x50','\xA4','\x40','\x57','\x1D','\x3E','\x1B',
    '\x05','\xA3','\x14','\xFC','\xC8','\x09','\xA2','\x6D','\x57','\xF0','\xD1','\x3D','\x9C','\xC6','\x14','\xCA',
    '\x18','\x9A','\x5C','\xB4','\x86','\xD5','\xF5','\x56','\xA3','\xD7','\x18','\x8A','\x2D','\x3E','\xFE','\x68',
    '\xE6','\x31','\xF1','\xE4','\xB9','\x26','\xD6','\x9B','\xED','\xAE','\xC3','\xF3','\xE5','\x2E','\x2E','\x4A',
    '\x23','\xAD','\x3C','\xD7','\x8C','\x81','\x35','\x15','\x19','\xF4','\xE6','\xFB','\x84','\xFD','\x13','\x63',
    '\xF6','\x9C','\x7C','\xF2','\x10','\xA2','\xDB','\xCB','\xC3','\xF1','\xD4','\xBC','\x01','\xDB','\x85','\xA3',
    '\xD3','\xC3','\x00','\x00','\x00'
};

TEST_F(alignment_file_input_f, decompression_by_filename_gz)
{
    seqan3::test::tmp_filename filename{"alignment_file_output_test.sam.gz"};

    {
        std::ofstream of{filename.get_path(), std::ios::binary};

        std::copy(input_gz.begin(), input_gz.end(), std::ostreambuf_iterator<char>{of});
    }

    seqan3::alignment_file_input fin{filename.get_path()};

    decompression_impl(*this, fin);
}

TEST_F(alignment_file_input_f, decompression_by_stream_gz)
{
    seqan3::alignment_file_input fin{std::istringstream{input_gz}, seqan3::format_sam{}};

    decompression_impl(*this, fin);
}

TEST_F(alignment_file_input_f, read_empty_gz_file)
{
    std::string empty_zipped_file
    {
        '\x1f', '\x8b', '\x08', '\x08', '\x5a', '\x07', '\x98', '\x5c',
        '\x00', '\x03', '\x66', '\x6f', '\x6f', '\x00', '\x03', '\x00',
        '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00'
    };
    seqan3::alignment_file_input fin{std::istringstream{empty_zipped_file}, seqan3::format_sam{}};

    EXPECT_TRUE(fin.begin() == fin.end());
}

std::string input_bgzf
{
    '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43',
    '\x02', '\x00', '\xF2', '\x00', '\x6D', '\x8E', '\x5D', '\x4B', '\xC3', '\x40', '\x10', '\x45', '\x9F', '\x6F',
    '\x7E', '\x45', '\x43', '\xC5', '\x6A', '\x8D', '\x31', '\xFB', '\x61', '\x0A', '\xF3', '\x94', '\x35', '\x81',
    '\x55', '\x30', '\x9B', '\xCA', '\x06', '\x5F', '\x25', '\xD4', '\xA8', '\x45', '\xB3', '\x2B', '\x51', '\x51',
    '\xFF', '\xBD', '\xE9', '\xEA', '\x8B', '\x50', '\x18', '\xB8', '\xDC', '\x61', '\xCE', '\x61', '\x8A', '\xCB',
    '\x0A', '\xB7', '\x86', '\x58', '\x9A', '\xC3', '\x36', '\xF4', '\xE1', '\x9E', '\x9D', '\xFF', '\x74', '\xD0',
    '\x0D', '\x39', '\xEF', '\xFA', '\xA8', '\xB0', '\x37', '\xB0', '\x86', '\xC6', '\xFE', '\x01', '\xD7', '\x86',
    '\x84', '\x8C', '\x8A', '\xB5', '\xC6', '\x55', '\x45', '\xAF', '\xA3', '\x7F', '\x64', '\x58', '\x1B', '\xDA',
    '\x78', '\xFF', '\x72', '\xB7', '\x6B', '\x63', '\x37', '\x44', '\x45', '\xD9', '\xA0', '\x7D', '\xDA', '\xBE',
    '\xCD', '\xA6', '\xE9', '\x66', '\x1B', '\x3F', '\x0C', '\xBD', '\x7B', '\x4F', '\xA3', '\xB1', '\xEF', '\xEE',
    '\x19', '\x24', '\xC3', '\x4E', '\xC3', '\x90', '\x33', '\x30', '\xCB', '\x6A', '\x56', '\xF1', '\xFA', '\x77',
    '\x93', '\x41', '\x64', '\x19', '\x54', '\xA9', '\x5B', '\xC4', '\xF3', '\xF9', '\x01', '\x94', '\xA5', '\x2D',
    '\x71', '\x98', '\x7A', '\x8A', '\x55', '\xA0', '\x39', '\x24', '\x0F', '\xB7', '\x1C', '\x39', '\xC7', '\x6A',
    '\x62', '\x27', '\xDE', '\xFE', '\xA3', '\xB5', '\x2E', '\x5B', '\x6D', '\x94', '\x0E', '\x8A', '\xC3', '\xC5',
    '\xD1', '\xF1', '\x12', '\x5F', '\xDF', '\x74', '\x41', '\x36', '\x11', '\x89', '\x4C', '\xCE', '\x83', '\x46',
    '\x40', '\x8A', '\x00', '\x09', '\xE4', '\xE2', '\xEF', '\x09', '\xB9', '\x47', '\xA6', '\xB5', '\xD2', '\xAD',
    '\x6A', '\x15', '\xE2', '\x78', '\x79', '\x92', '\x9C', '\xA6', '\x67', '\xD1', '\x0F', '\x72', '\x0E', '\xE3',
    '\xE8', '\x26', '\x01', '\x00', '\x00', '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00',
    '\xFF', '\x06', '\x00', '\x42', '\x43', '\x02', '\x00', '\x1B', '\x00', '\x03', '\x00', '\x00', '\x00', '\x00',
    '\x00', '\x00', '\x00', '\x00', '\x00'
};

TEST_F(alignment_file_input_f, decompression_by_filename_bgzf)
{
    seqan3::test::tmp_filename filename{"alignment_file_output_test.sam.bgzf"};

    {
        std::ofstream of{filename.get_path(), std::ios::binary};

        std::copy(input_bgzf.begin(), input_bgzf.end(), std::ostreambuf_iterator<char>{of});
    }

    seqan3::alignment_file_input fin{filename.get_path()};

    decompression_impl(*this, fin);
}

TEST_F(alignment_file_input_f, decompression_by_stream_bgzf)
{
    seqan3::alignment_file_input fin{std::istringstream{input_bgzf}, seqan3::format_sam{}};

    decompression_impl(*this, fin);
}

TEST_F(alignment_file_input_f, read_empty_bgzf_file)
{
    std::string empty_bgzf_file
    {
        '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF',
        '\x06', '\x00', '\x42', '\x43', '\x02', '\x00', '\x1B', '\x00', '\x03', '\x00',
        '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00'
    };
    seqan3::alignment_file_input fin{std::istringstream{empty_bgzf_file}, seqan3::format_sam{}};

    EXPECT_TRUE(fin.begin() == fin.end());
}

#endif

#ifdef SEQAN3_HAS_BZIP2
std::string input_bz2
{
    '\x42','\x5A','\x68','\x39','\x31','\x41','\x59','\x26','\x53','\x59','\x7B','\xE2','\xE1','\x92','\x00','\x00',
    '\x5C','\x5F','\x80','\x00','\x30','\x2D','\xFF','\xFF','\x90','\x3C','\x83','\x0C','\x00','\x27','\x20','\x10',
    '\x60','\x20','\x00','\x8A','\x86','\x82','\x4D','\x4D','\xA6','\xA6','\x9A','\x60','\xD4','\xC8','\xC9','\x99',
    '\x35','\x34','\x06','\x44','\x9B','\x51','\xA0','\x83','\x4C','\x06','\x88','\x31','\x19','\x32','\xDF','\x59',
    '\x81','\x84','\x10','\x62','\x4B','\x06','\x22','\x21','\xA8','\xEA','\x68','\xCD','\xA2','\x15','\xB7','\xE5',
    '\xA7','\xAB','\x0A','\xD2','\xB8','\x0A','\xEF','\xC3','\x18','\x35','\xFE','\x2C','\xE9','\x1C','\x72','\x8D',
    '\xA6','\xE2','\xC7','\x3D','\xBC','\x41','\x0E','\x00','\x50','\x3E','\x05','\x0E','\x0F','\x46','\xF5','\x2B',
    '\x39','\xEF','\x92','\xA5','\x28','\x85','\xEA','\xA5','\x93','\xE0','\xFD','\x27','\xBF','\x76','\xCC','\xE2',
    '\x6A','\xE9','\x32','\xE0','\x11','\x05','\x09','\x44','\xD2','\x51','\xB6','\x90','\x2A','\x73','\x94','\x54',
    '\x62','\x96','\x19','\xBB','\x92','\xE4','\xB8','\x20','\x28','\x32','\x8E','\x0C','\x09','\xE7','\xF8','\xBB',
    '\x92','\x29','\xC2','\x84','\x83','\xDF','\x17','\x0C','\x90'
};

TEST_F(alignment_file_input_f, decompression_by_filename_bz2)
{
    seqan3::test::tmp_filename filename{"alignment_file_output_test.sam.bz2"};

    {
        std::ofstream of{filename.get_path(), std::ios::binary};

        std::copy(input_bz2.begin(), input_bz2.end(), std::ostreambuf_iterator<char>{of});
    }

    seqan3::alignment_file_input fin{filename.get_path()};

    decompression_impl(*this, fin);
}

TEST_F(alignment_file_input_f, decompression_by_stream_bz2)
{
    seqan3::alignment_file_input fin{std::istringstream{input_bz2}, seqan3::format_sam{}};

    decompression_impl(*this, fin);
}

TEST_F(alignment_file_input_f, read_empty_bz2_file)
{
    std::string empty_zipped_file
    {
        '\x42', '\x5a', '\x68', '\x39', '\x17', '\x72', '\x45', '\x38', '\x50', '\x90', '\x00', '\x00', '\x00', '\x00'
    };
    seqan3::alignment_file_input fin{std::istringstream{empty_zipped_file}, seqan3::format_sam{}};

    EXPECT_TRUE(fin.begin() == fin.end());
}
#endif

// ----------------------------------------------------------------------------
// SAM format specificities
// ----------------------------------------------------------------------------

struct alignment_file_input_sam_format_f : public alignment_file_input_f
{
    std::vector<seqan3::dna4_vector> const ref_seqs = {"ACTGATCGAGAGGATCTAGAGGAGATCGTAGGAC"_dna4};
    std::vector<std::string> const ref_ids = {"ref"};

    std::vector<seqan3::gapped<seqan3::dna4>> ref_seq_gapped1 = {'A'_dna4, 'C'_dna4, 'T'_dna4, 'G'_dna4};
    std::vector<seqan3::gapped<seqan3::dna4>> ref_seq_gapped2 = {'C'_dna4, 'T'_dna4, 'G'_dna4, 'A'_dna4,
                                                                 'T'_dna4, 'C'_dna4, 'G'_dna4, 'A'_dna4,
                                                                 'G'_dna4};
    std::vector<seqan3::gapped<seqan3::dna4>> ref_seq_gapped3 = {'T'_dna4, 'G'_dna4, 'A'_dna4, 'T'_dna4,
                                                                 'C'_dna4, 'G'_dna4, 'A'_dna4, 'G'_dna4,};

    std::vector<std::pair<std::vector<seqan3::gapped<seqan3::dna4>>,
                          std::vector<seqan3::gapped<seqan3::dna5>>>> alignments_expected
    {
        {ref_seq_gapped1, std::vector<seqan3::gapped<seqan3::dna5>>{'C'_dna5, seqan3::gap{}, 'G'_dna5, 'T'_dna5}},
        {ref_seq_gapped2, std::vector<seqan3::gapped<seqan3::dna5>>{'A'_dna5, 'G'_dna5, 'G'_dna5, 'C'_dna5, 'T'_dna5,
                                                                    'G'_dna5, 'N'_dna5, seqan3::gap{}, 'A'_dna5}},
        {ref_seq_gapped3, std::vector<seqan3::gapped<seqan3::dna5>>{'G'_dna5, seqan3::gap{}, 'A'_dna5, 'G'_dna5,
                                                                    'T'_dna5, 'A'_dna5, seqan3::gap{}, 'T'_dna5}}
    };
};

TEST_F(alignment_file_input_sam_format_f, construct_by_filename_and_read_alignments)
{
    seqan3::test::tmp_filename filename{"alignment_file_input_constructor.sam"};
    {
        std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        filecreator << input;
    }

    seqan3::alignment_file_input fin{filename.get_path(), ref_ids, ref_seqs, seqan3::fields<seqan3::field::alignment>{}};

    EXPECT_EQ(fin.header().ref_ids(), ref_ids);

    size_t counter = 0;
    for (auto & [ alignment ] : fin)
    {
        EXPECT_TRUE(std::ranges::equal(std::get<0>(alignment), std::get<0>(alignments_expected[counter])));
        EXPECT_TRUE(std::ranges::equal(std::get<1>(alignment), std::get<1>(alignments_expected[counter])));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

TEST_F(alignment_file_input_sam_format_f, construct_from_stream_and_read_alignments)
{
    seqan3::alignment_file_input fin{std::istringstream{input},
                                     ref_ids,
                                     ref_seqs,
                                     seqan3::format_sam{},
                                     seqan3::fields<seqan3::field::alignment>{}};

    EXPECT_EQ(fin.header().ref_ids(), ref_ids);

    size_t counter = 0;
    for (auto & [ alignment ] : fin)
    {
        EXPECT_TRUE(std::ranges::equal(std::get<0>(alignment), std::get<0>(alignments_expected[counter])));
        EXPECT_TRUE(std::ranges::equal(std::get<1>(alignment), std::get<1>(alignments_expected[counter])));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

TEST_F(alignment_file_input_sam_format_f, construct_from_stream_and_read_alignment_with_dummy)
{
    seqan3::alignment_file_input fin{std::istringstream{input},
                                     seqan3::format_sam{},
                                     seqan3::fields<seqan3::field::alignment>{}};

    size_t counter = 0;
    for (auto & [ alignment ] : fin)
    {
        EXPECT_EQ(std::get<1>(alignment), std::get<1>(alignments_expected[counter]));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

// ----------------------------------------------------------------------------
// BAM format specificities
// ----------------------------------------------------------------------------

struct alignment_file_input_bam_format_f : public alignment_file_input_sam_format_f
{
    std::string binary_input{ // corresponds to 'input' from alignment_file_input_f fixture
        '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42',
        '\x43', '\x02', '\x00', '\x8D', '\x00', '\x73', '\x72', '\xF4', '\x65', '\x4C', '\x66', '\x60', '\x60',
        '\x70', '\xF0', '\x70', '\xE1', '\x0C', '\xF3', '\xB3', '\x32', '\xD4', '\x33', '\xE3', '\x0C', '\xF6',
        '\xB7', '\x2A', '\xCD', '\xCB', '\xCE', '\xCB', '\x2F', '\xCF', '\xE3', '\x74', '\xF7', '\xB7', '\xCA',
        '\xCB', '\xCF', '\x4B', '\xE5', '\x72', '\x08', '\x0E', '\xE4', '\x0C', '\xF6', '\xB3', '\x2A', '\x4A',
        '\x4D', '\xE3', '\xF4', '\xF1', '\xB3', '\x32', '\x36', '\xE1', '\x72', '\x08', '\x70', '\xE7', '\xF4',
        '\x74', '\xB1', '\x2A', '\x28', '\xCA', '\x4F', '\x37', '\xE4', '\x0C', '\xF0', '\xB3', '\x4A', '\xCE',
        '\xCF', '\xCF', '\x89', '\x07', '\xF1', '\x8A', '\x12', '\x73', '\xB9', '\x1C', '\x9C', '\xFD', '\x39',
        '\x43', '\x32', '\x32', '\x8B', '\x15', '\x80', '\x28', '\x51', '\x21', '\x39', '\x3F', '\x37', '\x37',
        '\x35', '\xAF', '\x44', '\x8F', '\x8B', '\x11', '\x68', '\x0D', '\x0B', '\x10', '\x03', '\x4D', '\x61',
        '\x50', '\x02', '\xD2', '\x00', '\xF2', '\x5C', '\x8E', '\x8D', '\x7B', '\x00', '\x00', '\x00', '\x1F',
        '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43',
        '\x02', '\x00', '\x9F', '\x00', '\x73', '\x61', '\x40', '\x00', '\x36', '\x5B', '\x4F', '\x21', '\x16',
        '\x06', '\x4D', '\x06', '\x16', '\x28', '\x9F', '\x13', '\x88', '\x75', '\x18', '\x19', '\x18', '\x8A',
        '\x52', '\x13', '\x53', '\x0C', '\x19', '\x44', '\x80', '\x3C', '\x01', '\x20', '\x16', '\x02', '\x62',
        '\x05', '\x10', '\xED', '\xC1', '\xC0', '\xC4', '\xC4', '\xEC', '\x18', '\xEC', '\xCC', '\xE4', '\xE7',
        '\xEB', '\xCC', '\x1E', '\x04', '\xD5', '\xC3', '\x08', '\x32', '\xC7', '\x0E', '\x64', '\x8E', '\x16',
        '\x58', '\x3F', '\xBA', '\x39', '\x46', '\x0C', '\x05', '\x50', '\x33', '\x40', '\x66', '\x81', '\xCC',
        '\x14', '\x71', '\x6A', '\xF9', '\xE8', '\x00', '\x32', '\x8A', '\x95', '\x8D', '\x9D', '\x83', '\xB3',
        '\xA2', '\xD2', '\x29', '\x98', '\x19', '\x28', '\xCA', '\x0C', '\x74', '\x05', '\x2B', '\x83', '\x1F',
        '\xD4', '\x04', '\x26', '\x90', '\xA9', '\xF6', '\x9E', '\x42', '\xEC', '\x0C', '\xDA', '\x0C', '\x1C',
        '\x58', '\x4C', '\x35', '\x46', '\x71', '\x9D', '\x03', '\x9A', '\x0D', '\x2E', '\x22', '\x8D', '\x8D',
        '\x40', '\xF5', '\x5C', '\xDC', '\x3C', '\xBC', '\x7C', '\x00', '\xC5', '\xFD', '\x4B', '\xCD', '\xF0',
        '\x00', '\x00', '\x00', '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF',
        '\x06', '\x00', '\x42', '\x43', '\x02', '\x00', '\x1B', '\x00', '\x03', '\x00', '\x00', '\x00', '\x00',
        '\x00', '\x00', '\x00', '\x00', '\x00'
    };
};

#if SEQAN3_HAS_ZLIB
TEST_F(alignment_file_input_bam_format_f, construct_by_filename)
{
    seqan3::test::tmp_filename filename{"alignment_file_input_constructor.bam"};
    {
        std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        filecreator << binary_input;
    }

    seqan3::alignment_file_input fin{filename.get_path(), ref_ids, ref_seqs, seqan3::fields<seqan3::field::id,
                                                                                            seqan3::field::seq,
                                                                                            seqan3::field::qual,
                                                                                            seqan3::field::alignment>{}};

    EXPECT_EQ(fin.header().ref_ids(), ref_ids);
    EXPECT_EQ(fin.header().comments[0], std::string{"This is a comment."});

    size_t counter = 0;
    for (auto & [ id, seq, qual, alignment ] : fin)
    {
        EXPECT_EQ(id, id_comp[counter]);
        EXPECT_EQ(seq, seq_comp[counter]);
        EXPECT_EQ(qual, qual_comp[counter]);

        EXPECT_TRUE(std::ranges::equal(std::get<0>(alignment), std::get<0>(alignments_expected[counter])));
        EXPECT_TRUE(std::ranges::equal(std::get<1>(alignment), std::get<1>(alignments_expected[counter])));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

TEST_F(alignment_file_input_bam_format_f, construct_by_stream)
{
    std::istringstream stream{binary_input};

    seqan3::alignment_file_input fin{stream,
                                     ref_ids,
                                     ref_seqs,
                                     seqan3::format_bam{},
                                     seqan3::fields<seqan3::field::id,
                                                    seqan3::field::seq,
                                                    seqan3::field::qual,
                                                    seqan3::field::alignment>{}};

    EXPECT_EQ(fin.header().ref_ids(), ref_ids);
    EXPECT_EQ(fin.header().comments[0], std::string{"This is a comment."});

    size_t counter = 0;
    for (auto & [ id, seq, qual, alignment ] : fin)
    {
        EXPECT_EQ(id, id_comp[counter]);
        EXPECT_EQ(seq, seq_comp[counter]);
        EXPECT_EQ(qual, qual_comp[counter]);

        EXPECT_TRUE(std::ranges::equal(std::get<0>(alignment), std::get<0>(alignments_expected[counter])));
        EXPECT_TRUE(std::ranges::equal(std::get<1>(alignment), std::get<1>(alignments_expected[counter])));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}
#endif // SEQAN3_HAS_ZLIB
