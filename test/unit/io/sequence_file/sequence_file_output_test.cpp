// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iterator>
#include <sstream>

#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/test/tmp_directory.hpp>
#include <seqan3/test/zlib_skip.hpp>
#include <seqan3/utility/views/zip.hpp>

using seqan3::operator""_dna5;
using seqan3::operator""_phred42;

using default_fields = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::qual>;

TEST(sequence_file_output_iterator, concepts)
{
    using it_t = typename seqan3::sequence_file_output<>::iterator;
    using sen_t = typename seqan3::sequence_file_output<>::sentinel;

    EXPECT_TRUE((std::output_iterator<it_t, std::tuple<std::string, std::string>>));
    EXPECT_TRUE((std::sentinel_for<sen_t, it_t>));
}

std::vector<seqan3::dna5_vector> seqs{
    "ACGT"_dna5,
    "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN"_dna5,
    "GGAGTATAATATATATATATATAT"_dna5};

std::vector<std::string> ids{"TEST 1", "Test2", "Test3"};

std::string const output_comp{
    ">TEST 1\n"
    "ACGT\n"
    ">Test2\n"
    "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN\n"
    ">Test3\n"
    "GGAGTATAATATATATATATATAT\n"};

std::vector<std::vector<seqan3::phred42>> quals{
    "!!!!"_phred42,
    "!#@$!#@$!#@#!$@#!$@#!$!#@$!#@#!$@#!$!#$@!!$$$$$$$$$$$$!!!!!!!!!!!!!!!!!!!!$$$$$$$$$$!!!!!$!"_phred42,
    "!@#!@#!#!######@$!#@!!!@"_phred42};

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TEST(general, concepts)
{
    using t = seqan3::sequence_file_output<>;
    EXPECT_TRUE((std::ranges::output_range<t, std::tuple<std::string, std::string>>));

    using ct = seqan3::sequence_file_output<> const;
    // not const-iterable
    EXPECT_FALSE((std::ranges::output_range<ct, std::tuple<std::string, std::string>>));
}

TEST(general, construct_by_filename)
{
    /* just the filename */
    {
        seqan3::test::tmp_directory tmp;
        auto filename = tmp.path() / "sequence_file_output_constructor.fasta";
        EXPECT_NO_THROW(seqan3::sequence_file_output<>{filename});
    }

    /* wrong extension */
    {
        seqan3::test::tmp_directory tmp;
        auto filename = tmp.path() / "sequence_file_output_constructor.xyz";
        std::ofstream filecreator{filename, std::ios::out | std::ios::binary};
        EXPECT_THROW(seqan3::sequence_file_output<>{filename}, seqan3::unhandled_extension_error);
    }

    /* unknown file */
    {
        seqan3::test::tmp_directory tmp;
        auto filename = tmp.path() / "I/do/not/exist.fasta";
        EXPECT_THROW(seqan3::sequence_file_output<>{filename}, seqan3::file_open_error);
    }

    /* filename + fields */
    using fields_seq = seqan3::fields<seqan3::field::seq>;
    {
        seqan3::test::tmp_directory tmp;
        auto filename = tmp.path() / "sequence_file_output_constructor.fasta";
        EXPECT_NO_THROW((
            seqan3::sequence_file_output<fields_seq, seqan3::type_list<seqan3::format_fasta>>{filename, fields_seq{}}));
    }
}

TEST(general, construct_from_stream)
{
    /* stream + format_tag */
    EXPECT_NO_THROW((
        seqan3::sequence_file_output<default_fields, seqan3::type_list<seqan3::format_fasta>>{std::ostringstream{},
                                                                                              seqan3::format_fasta{}}));

    /* stream + format_tag + fields */
    EXPECT_NO_THROW(
        (seqan3::sequence_file_output<default_fields, seqan3::type_list<seqan3::format_fasta>>{std::ostringstream{},
                                                                                               seqan3::format_fasta{},
                                                                                               default_fields{}}));
}

TEST(general, default_template_args_and_deduction_guides)
{
    using comp2 = seqan3::type_list<seqan3::format_embl,
                                    seqan3::format_fasta,
                                    seqan3::format_fastq,
                                    seqan3::format_genbank,
                                    seqan3::format_sam>;
    using comp3 = char;

    /* default template args */
    {
        using t = seqan3::sequence_file_output<>;
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, default_fields>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats, comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
    }

    /* guided filename constructor */
    {
        seqan3::test::tmp_directory tmp;
        auto filename = tmp.path() / "sequence_file_output_constructor.fasta";

        seqan3::sequence_file_output fout{filename};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, default_fields>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats, comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
    }

    /* guided filename constructor + custom fields */
    {
        seqan3::test::tmp_directory tmp;
        auto filename = tmp.path() / "sequence_file_output_constructor.fasta";

        seqan3::sequence_file_output fout{filename, seqan3::fields<seqan3::field::seq>{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, seqan3::fields<seqan3::field::seq>>)); // changed
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats, comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
    }

    /* guided stream constructor */
    {
        std::ostringstream ext{};
        seqan3::sequence_file_output fout{ext, seqan3::format_fasta{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, default_fields>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats, seqan3::type_list<seqan3::format_fasta>>)); // changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
    }

    /* guided stream temporary constructor */
    {
        seqan3::sequence_file_output fout{std::ostringstream{}, seqan3::format_fasta{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, default_fields>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats, seqan3::type_list<seqan3::format_fasta>>)); // changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_char_type, comp3>));
    }
}

// ----------------------------------------------------------------------------
// *impl
// ----------------------------------------------------------------------------

template <typename fn_t>
void row_wise_impl(fn_t fn)
{
    seqan3::sequence_file_output fout{std::ostringstream{}, seqan3::format_fasta{}};
    fout.options.fasta_letters_per_line = 0;

    for (size_t i = 0; i < 3; ++i)
        fn(fout, i);

    fout.get_stream().flush();
    EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), output_comp);
}

template <typename source_t>
void assign_impl(source_t && source)
{
    seqan3::sequence_file_output fout{std::ostringstream{}, seqan3::format_fasta{}};
    fout.options.fasta_letters_per_line = 0;

    fout = source;

    fout.get_stream().flush();
    EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), output_comp);
}

// ----------------------------------------------------------------------------
// row
// ----------------------------------------------------------------------------

TEST(row, assign_to_iterator)
{
    using fields_seq_id = seqan3::fields<seqan3::field::seq, seqan3::field::id>;

    row_wise_impl(
        [&](auto & file, size_t i)
        {
            seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>, fields_seq_id> r{seqs[i], ids[i]};

            std::ranges::begin(file) = r;
        });
}

TEST(row, push_back_record)
{
    using fields_seq_id = seqan3::fields<seqan3::field::seq, seqan3::field::id>;

    row_wise_impl(
        [&](auto & file, size_t i)
        {
            seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>, fields_seq_id> r{seqs[i], ids[i]};

            file.push_back(r);
        });
}

TEST(row, push_back_record_rvalue)
{
    using fields_seq_id = seqan3::fields<seqan3::field::seq, seqan3::field::id>;

    row_wise_impl(
        [&](auto & file, size_t i)
        {
            seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>, fields_seq_id> r{seqs[i], ids[i]};

            file.push_back(std::move(r));
        });
}

TEST(row, push_back_record_const)
{
    using fields_seq_id = seqan3::fields<seqan3::field::seq, seqan3::field::id>;

    row_wise_impl(
        [&](auto & file, size_t i)
        {
            seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>, fields_seq_id> const r{seqs[i], ids[i]};

            file.push_back(r);
        });
}

TEST(row, push_back_record_const_element)
{
    using fields_seq_id = seqan3::fields<seqan3::field::seq, seqan3::field::id>;

    row_wise_impl(
        [&](auto & file, size_t i)
        {
            seqan3::record<seqan3::type_list<seqan3::dna5_vector const, std::string const>, fields_seq_id> const r{
                seqs[i],
                ids[i]};

            file.push_back(r);
        });
}

TEST(row, push_back_tuple)
{
    row_wise_impl(
        [&](auto & file, size_t i)
        {
            std::tuple<seqan3::dna5_vector, std::string> t{seqs[i], ids[i]};

            file.push_back(t);
        });
}

TEST(row, push_back_tuple_rvalue)
{
    row_wise_impl(
        [&](auto & file, size_t i)
        {
            std::tuple<seqan3::dna5_vector, std::string> t{seqs[i], ids[i]};

            file.push_back(std::move(t));
        });
}

TEST(row, push_back_tuple_const)
{
    row_wise_impl(
        [&](auto & file, size_t i)
        {
            std::tuple<seqan3::dna5_vector, std::string> const t{seqs[i], ids[i]};

            file.push_back(t);
        });
}

TEST(row, push_back_tuple_const_element)
{
    row_wise_impl(
        [&](auto & file, size_t i)
        {
            std::tuple<seqan3::dna5_vector const, std::string const> t{seqs[i], ids[i]};

            file.push_back(t);
        });
}

TEST(row, emplace_back)
{
    row_wise_impl(
        [&](auto & file, size_t i)
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
                   seqan3::fields<seqan3::field::qual, seqan3::field::id, seqan3::field::seq>>
        rec{qual, ids[1], seqs[1]};

    seqan3::sequence_file_output fout{std::ostringstream{},
                                      seqan3::format_fasta{},
                                      seqan3::fields<seqan3::field::seq, seqan3::field::id>{}};
    fout.push_back(rec);
    fout.get_stream().flush();

    std::string const expected_out{">Test2\n"
                                   "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGG\n"
                                   "CTGNAGGCTGN\n"};
    EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), expected_out);
}

// ----------------------------------------------------------------------------
// rows
// ----------------------------------------------------------------------------

TEST(rows, assign_range_of_records)
{
    std::vector<seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>,
                               seqan3::fields<seqan3::field::seq, seqan3::field::id>>>
        range;

    for (size_t i = 0; i < 3; ++i)
        range.emplace_back(seqs[i], ids[i]);

    assign_impl(range);
}

TEST(rows, assign_range_of_records_const)
{
    std::vector<seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>,
                               seqan3::fields<seqan3::field::seq, seqan3::field::id>>>
        range;

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

// ----------------------------------------------------------------------------
// columns
// ----------------------------------------------------------------------------

TEST(columns, assign_tuple_of_columns)
{
    assign_impl(seqan3::views::zip(seqs, ids));
}

TEST(columns, writing_id_seq_qual)
{
    seqan3::sequence_file_output fout{std::ostringstream{},
                                      seqan3::format_fasta{},
                                      seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::qual>()};
    fout.options.fasta_letters_per_line = 0;

    fout = seqan3::views::zip(ids, seqs, quals);

    fout.get_stream().flush();
    EXPECT_EQ(reinterpret_cast<std::ostringstream &>(fout.get_stream()).str(), output_comp);
}

// ----------------------------------------------------------------------------
// compression
// ----------------------------------------------------------------------------

std::string compression_by_filename_impl(seqan3::test::sandboxed_path const & filename)
{
    {
        seqan3::sequence_file_output fout{filename};
        fout.options.fasta_blank_before_id = true;
        fout.options.fasta_letters_per_line = 0;

        for (size_t i = 0; i < 3; ++i)
        {
            seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>,
                           seqan3::fields<seqan3::field::seq, seqan3::field::id>>
                r{seqs[i], ids[i]};

            fout.push_back(r);
        }
    }

    std::string buffer;

    {
        std::ifstream fi{filename, std::ios::binary};

        buffer = std::string{std::istreambuf_iterator<char>{fi}, std::istreambuf_iterator<char>{}};
    }
    return buffer;
}

template <typename comp_stream_t>
void compression_by_stream_impl(comp_stream_t & stream)
{
    seqan3::sequence_file_output fout{stream, seqan3::format_fasta{}};
    fout.options.fasta_blank_before_id = true;
    fout.options.fasta_letters_per_line = 0;

    for (size_t i = 0; i < 3; ++i)
    {
        seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>,
                       seqan3::fields<seqan3::field::seq, seqan3::field::id>>
            r{seqs[i], ids[i]};

        fout.push_back(r);
    }
}

#if defined(SEQAN3_HAS_ZLIB)
std::string expected_gz{'\x1F', '\x8B', '\x08', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\xB3', '\x53',
                        '\x08', '\x71', '\x0D', '\x0E', '\x51', '\x30', '\xE4', '\x72', '\x74', '\x76', '\x0F', '\xE1',
                        '\xB2', '\x53', '\x08', '\x49', '\x2D', '\x2E', '\x31', '\xE2', '\x72', '\x74', '\x77', '\x77',
                        '\x0E', '\x71', '\xF7', '\xA3', '\x05', '\x05', '\xB5', '\xC3', '\x98', '\xCB', '\xDD', '\xDD',
                        '\xD1', '\x3D', '\xC4', '\x31', '\xC4', '\xD1', '\x31', '\x04', '\x15', '\x72', '\x01', '\x00',
                        '\x27', '\xAD', '\xB4', '\xE9', '\x93', '\x00', '\x00', '\x00'};

std::string expected_bgzf{
    '\x1F', '\x8B', '\x08', '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x06', '\x00', '\x42',
    '\x43', '\x02', '\x00', '\x4A', '\x00', '\xB3', '\x53', '\x08', '\x71', '\x0D', '\x0E', '\x51', '\x30',
    '\xE4', '\x72', '\x74', '\x76', '\x0F', '\xE1', '\xB2', '\x53', '\x08', '\x49', '\x2D', '\x2E', '\x31',
    '\xE2', '\x72', '\x74', '\x77', '\x77', '\x0E', '\x71', '\xF7', '\xA3', '\x05', '\x05', '\xB5', '\xC3',
    '\x98', '\xCB', '\xDD', '\xDD', '\xD1', '\x3D', '\xC4', '\x31', '\xC4', '\x11', '\x88', '\x50', '\x20',
    '\x17', '\x00', '\x27', '\xAD', '\xB4', '\xE9', '\x93', '\x00', '\x00', '\x00', '\x1F', '\x8B', '\x08',
    '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xFF', '\x06', '\x00', '\x42', '\x43', '\x02', '\x00',
    '\x1B', '\x00', '\x03', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00', '\x00'};

TEST(compression, by_filename_gz)
{
    seqan3::test::tmp_directory tmp;
    auto filename = tmp.path() / "sequence_file_output_test.fasta.gz";

    std::string buffer = compression_by_filename_impl(filename);
    buffer[9] = '\x00'; // zero out OS byte
    SEQAN3_TEST_GTEST_SKIP_ZLIB_DEFLATE;
    EXPECT_EQ(buffer, expected_gz);
}

TEST(compression, by_stream_gz)
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

TEST(compression, by_filename_bgzf)
{
    seqan3::test::tmp_directory tmp;
    auto filename = tmp.path() / "sequence_file_output_test.fasta.bgzf";

    std::string buffer = compression_by_filename_impl(filename);
    buffer[9] = '\x00'; // zero out OS byte
    SEQAN3_TEST_GTEST_SKIP_ZLIB_DEFLATE;
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
    buffer[9] = '\x00'; // zero out OS byte
    SEQAN3_TEST_GTEST_SKIP_ZLIB_DEFLATE;
    EXPECT_EQ(buffer, expected_bgzf);
}

#endif

#if defined(SEQAN3_HAS_BZIP2)
std::string expected_bz2{'\x42', '\x5A', '\x68', '\x39', '\x31', '\x41', '\x59', '\x26', '\x53', '\x59', '\xB4',
                         '\x68', '\xEA', '\xE3', '\x00', '\x00', '\x06', '\xDF', '\x80', '\x00', '\x10', '\x40',
                         '\x00', '\x38', '\x01', '\x2A', '\x81', '\x0C', '\x00', '\x02', '\x00', '\x0C', '\x00',
                         '\x20', '\x00', '\x50', '\xA6', '\x00', '\x09', '\xA0', '\x8A', '\x10', '\x9A', '\x32',
                         '\x34', '\xD9', '\xAB', '\x5F', '\x16', '\xE9', '\xEB', '\x86', '\x5B', '\x46', '\x41',
                         '\x8D', '\xD0', '\x1E', '\x12', '\x8C', '\xC0', '\xB5', '\x48', '\xD2', '\x3A', '\x9B',
                         '\x23', '\xB9', '\x9F', '\x64', '\x98', '\x1E', '\xEE', '\x8C', '\x18', '\x3E', '\x38',
                         '\x7E', '\x2E', '\xE4', '\x8A', '\x70', '\xA1', '\x21', '\x68', '\xD1', '\xD5', '\xC6'};

TEST(compression, by_filename_bz2)
{
    seqan3::test::tmp_directory tmp;
    auto filename = tmp.path() / "sequence_file_output_test.fasta.bz2";

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
