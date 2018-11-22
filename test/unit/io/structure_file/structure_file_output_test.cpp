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

/*!\file
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \brief Tests for writing sequence files with structure information.
 */

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
#include <seqan3/std/view/filter.hpp>
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

    /* non-existent file*/
    {
        EXPECT_THROW(structure_file_out<>{"/dev/nonexistant/foobarOOO"}, file_open_error);
    }

    /* filename + fields */
    {
        test::tmp_filename filename{"structure_file_output_constructor.dbn"};
        EXPECT_NO_THROW((structure_file_out<fields<field::SEQ>,
                                            type_list<structure_file_format_vienna>,
                                            std::ofstream>
                                            {filename.get_path(), fields<field::SEQ>{}}));
    }
}

TEST(structure_file_output_class, construct_from_stream)
{
    /* stream + format_tag */
    EXPECT_NO_THROW((structure_file_out<fields<field::SEQ, field::ID, field::STRUCTURE>,
                                        type_list<structure_file_format_vienna>,
                                        std::ostringstream>
                                        {std::ostringstream{}, structure_file_format_vienna{}}));

    /* stream + format_tag + fields */
    EXPECT_NO_THROW((structure_file_out<fields<field::SEQ, field::ID, field::STRUCTURE>,
                                        type_list<structure_file_format_vienna>,
                                        std::ostringstream>
                     {std::ostringstream{}, structure_file_format_vienna{},
                      fields<field::SEQ, field::ID, field::STRUCTURE>{}}));
}

TEST(structure_file_output_class, default_template_args_and_deduction_guides)
{
    using comp1 = fields<field::SEQ, field::ID, field::STRUCTURE>;
    using comp2 = type_list<structure_file_format_vienna>;
    using comp3 = std::ofstream;

    /* default template args */
    {
        using t = structure_file_out<>;
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
    }

    /* guided filename constructor */
    {
        test::tmp_filename filename{"structure_file_output_constructor.dbn"};
        structure_file_out fout{filename.get_path()};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
    }

    /* guided filename constructor + custom fields */
    {
        test::tmp_filename filename{"structure_file_output_constructor.dbn"};
        structure_file_out fout{filename.get_path(), fields<field::SEQ>{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
    }

    /* guided stream constructor */
    {
        structure_file_out fout{std::ostringstream{}, structure_file_format_vienna{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<structure_file_format_vienna>>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        std::ostringstream>));
    }

    /* guided stream constructor + custom fields */
    {
        structure_file_out fout{std::ostringstream{}, structure_file_format_vienna{}, fields<field::SEQ>{}};

        using t = decltype(fout);
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<structure_file_format_vienna>>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        std::ostringstream>));
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
        EXPECT_EQ(fout.get_stream().str(), output_comp);
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
        EXPECT_EQ(fout.get_stream().str(), output_comp);
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
    EXPECT_EQ(fout.get_stream().str(), output_comp);
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
