// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/view/zip.hpp>
#include <range/v3/view/filter.hpp>

#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/range/view/convert.hpp>
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

    /* filename + fields */
    {
        test::tmp_filename filename{"sequence_file_output_constructor.fasta"};
        EXPECT_NO_THROW(( sequence_file_output<fields<field::SEQ>,
                                            type_list<sequence_file_format_fasta>,
                                            std::ofstream>{filename.get_path(), fields<field::SEQ>{}} ));
    }
}

TEST(general, construct_from_stream)
{
    /* stream + format_tag */
    EXPECT_NO_THROW(( sequence_file_output<fields<field::SEQ, field::ID, field::QUAL>,
                                        type_list<sequence_file_format_fasta>,
                                        std::ostringstream>{std::ostringstream{},
                                                            sequence_file_format_fasta{}} ));


    /* stream + format_tag + fields */
    EXPECT_NO_THROW(( sequence_file_output<fields<field::SEQ, field::ID, field::QUAL>,
                                        type_list<sequence_file_format_fasta>,
                                        std::ostringstream>{std::ostringstream{},
                                                            sequence_file_format_fasta{},
                                                            fields<field::SEQ, field::ID, field::QUAL>{}} ));
}

TEST(general, default_template_args_and_deduction_guides)
{
    using comp1 = fields<field::SEQ, field::ID, field::QUAL>;
    using comp2 = type_list<sequence_file_format_fasta, sequence_file_format_fastq>;
    using comp3 = std::ofstream;

    /* default template args */
    {
        using t = sequence_file_output<>;
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
    }

    /* guided filename constructor */
    {
        test::tmp_filename filename{"sequence_file_output_constructor.fasta"};

        sequence_file_output fout{filename.get_path()};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
    }

    /* guided filename constructor + custom fields */
    {
        test::tmp_filename filename{"sequence_file_output_constructor.fasta"};

        sequence_file_output fout{filename.get_path(), fields<field::SEQ>{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));                   // changed
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
    }

    /* guided stream constructor */
    {
        sequence_file_output fout{std::ostringstream{}, sequence_file_format_fasta{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<sequence_file_format_fasta,
                                                                              sequence_file_format_fastq>>));// changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        std::ostringstream>));                   // changed
    }

    /* guided stream constructor + custom fields */
    {
        sequence_file_output fout{std::ostringstream{}, sequence_file_format_fasta{}, fields<field::SEQ>{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));                   // changed
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<sequence_file_format_fasta>>));// changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        std::ostringstream>));                   // changed
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
    EXPECT_EQ(fout.get_stream().str(), output_comp);
}

template <typename source_t>
void assign_impl(source_t && source)
{
    sequence_file_output fout{std::ostringstream{}, sequence_file_format_fasta{}};
    fout.options.fasta_letters_per_line = 0;

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

    sequence_file_output fout{std::ostringstream{}, sequence_file_format_fasta{}, fields<field::SEQ, field::ID>{}};
    fout.push_back(rec);
    fout.get_stream().flush();

    std::string const expected_out
    {
        "> Test2\n"
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGG\n"
        "CTGNAGGCTGN\n"
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
