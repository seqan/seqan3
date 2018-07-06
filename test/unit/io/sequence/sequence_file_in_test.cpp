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

#include <seqan3/io/sequence/sequence_file_in.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/concept/iterator.hpp>
#include <seqan3/std/view/filter.hpp>
#include <seqan3/test/tmp_filename.hpp>

using namespace seqan3;
using namespace seqan3::literal;

TEST(sequence_file_in_iterator, concepts)
{
    using it_t = typename sequence_file_in<>::iterator;
    using sen_t = typename sequence_file_in<>::sentinel;

    EXPECT_TRUE((seqan3::input_iterator_concept<it_t>));
    EXPECT_TRUE((seqan3::sentinel_concept<sen_t, it_t>));
}

struct sequence_file_in_f : public ::testing::Test
{
    std::string input
    {
        "> TEST1\n"
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
        "TEST1",
        "Test2",
        "Test3"
    };
};

TEST_F(sequence_file_in_f, concepts)
{
    using t = sequence_file_in<>;
    EXPECT_TRUE((seqan3::input_range_concept<t>));

    using ct = sequence_file_in<> const;
    // not const-iterable
    EXPECT_FALSE((seqan3::input_range_concept<ct>));
}

TEST_F(sequence_file_in_f, construct_by_filename)
{
    /* just the filename */
    {
        test::tmp_filename filename{"sequence_file_in_constructor.fasta"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
            filecreator << "> ID\nACGT\n"; // must contain at least one record
        }

        EXPECT_NO_THROW( sequence_file_in<>{filename.get_path()} );
    }

    // correct format check is done by tests of that format

    /* wrong extension */
    {
        test::tmp_filename filename{"sequence_file_in_constructor.xyz"};
        std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
        EXPECT_THROW( sequence_file_in<>{filename.get_path()} ,
                      unhandled_extension_error );
    }

    /* non-existent file*/
    {
        EXPECT_THROW( sequence_file_in<>{"/dev/nonexistant/foobarOOO"},
                      file_open_error);
    }

    /* filename + fields */
    {
        test::tmp_filename filename{"sequence_file_in_constructor.fasta"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
            filecreator << "> ID\nACGT\n"; // must contain at least one record
        }

        EXPECT_NO_THROW(( sequence_file_in<sequence_file_in_default_traits_dna,
                                           fields<field::SEQ>,
                                           type_list<sequence_file_format_fasta>,
                                           std::ifstream>{filename.get_path(), fields<field::SEQ>{}} ));
    }
}

TEST_F(sequence_file_in_f, construct_from_stream)
{
    /* stream + format_tag */
    EXPECT_NO_THROW(( sequence_file_in<sequence_file_in_default_traits_dna,
                                       fields<field::SEQ, field::ID, field::QUAL>,
                                       type_list<sequence_file_format_fasta>,
                                       std::istringstream>{std::istringstream{input},
                                                           sequence_file_format_fasta{}} ));


    /* stream + format_tag + fields */
    EXPECT_NO_THROW(( sequence_file_in<sequence_file_in_default_traits_dna,
                                       fields<field::SEQ, field::ID, field::QUAL>,
                                       type_list<sequence_file_format_fasta>,
                                       std::istringstream>{std::istringstream{input},
                                                           sequence_file_format_fasta{},
                                                           fields<field::SEQ, field::ID, field::QUAL>{}} ));
}

TEST_F(sequence_file_in_f, default_template_args_and_deduction_guides)
{
    using comp0 = sequence_file_in_default_traits_dna;
    using comp1 = fields<field::SEQ, field::ID, field::QUAL>;
    using comp2 = type_list<sequence_file_format_fasta>;
    using comp3 = std::ifstream;

    /* default template args */
    {
        using t = sequence_file_in<>;
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
    }

    /* guided filename constructor */
    {
        test::tmp_filename filename{"sequence_file_in_constructor.fasta"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
            filecreator << "> ID\nACGT\n"; // must contain at least one record
        }

        sequence_file_in fin{filename.get_path()};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
    }

    /* guided filename constructor + custom fields */
    {
        test::tmp_filename filename{"sequence_file_in_constructor.fasta"};

        {
            std::ofstream filecreator{filename.get_path(), std::ios::out | std::ios::binary};
            filecreator << "> ID\nACGT\n"; // must contain at least one record
        }

        sequence_file_in fin{filename.get_path(), fields<field::SEQ>{}};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));                   // changed
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      comp2>));
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        comp3>));
    }

    /* guided stream constructor */
    {
        sequence_file_in fin{std::istringstream{input}, sequence_file_format_fasta{}};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, comp1>));
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<sequence_file_format_fasta>>));// changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        std::istringstream>));                   // changed
    }

    /* guided stream constructor + custom fields */
    {
        sequence_file_in fin{std::istringstream{input}, sequence_file_format_fasta{}, fields<field::SEQ>{}};

        using t = decltype(fin);
        EXPECT_TRUE((std::is_same_v<typename t::traits_type,        comp0>));
        EXPECT_TRUE((std::is_same_v<typename t::selected_field_ids, fields<field::SEQ>>));                   // changed
        EXPECT_TRUE((std::is_same_v<typename t::valid_formats,      type_list<sequence_file_format_fasta>>));// changed
        EXPECT_TRUE((std::is_same_v<typename t::stream_type,        std::istringstream>));                   // changed
    }
}

TEST_F(sequence_file_in_f, record_reading)
{
    /* record based reading */
    sequence_file_in fin{std::istringstream{input}, sequence_file_format_fasta{}};

    size_t counter = 0;
    for (auto & rec : fin)
    {
        std::cout << (get<field::SEQ>(rec) | view::to_char) << std::endl;
        std::cout << (get<field::ID>(rec) | view::to_char) << std::endl;
        EXPECT_TRUE((ranges::equal(get<field::SEQ>(rec), seq_comp[counter])));
        EXPECT_TRUE((ranges::equal(get<field::ID>(rec),  id_comp[counter])));
        EXPECT_TRUE(empty(get<field::QUAL>(rec)));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

TEST_F(sequence_file_in_f, record_reading_struct_bind)
{
    /* record based reading */
    sequence_file_in fin{std::istringstream{input}, sequence_file_format_fasta{}};

    size_t counter = 0;
    for (auto & [ seq, id, qual ] : fin)
    {
        EXPECT_TRUE((ranges::equal(seq, seq_comp[counter])));
        EXPECT_TRUE((ranges::equal(id,  id_comp[counter])));
        EXPECT_TRUE(empty(qual));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

TEST_F(sequence_file_in_f, record_reading_custom_fields)
{
    /* record based reading */
    sequence_file_in fin{std::istringstream{input},
                         sequence_file_format_fasta{},
                         fields<field::ID, field::SEQ_QUAL>{}};

    size_t counter = 0;
    for (auto & [ id, seq_qual ] : fin)
    {
        EXPECT_TRUE((ranges::equal(seq_qual | view::convert<dna5>, seq_comp[counter])));
        EXPECT_TRUE((ranges::equal(id,  id_comp[counter])));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

TEST_F(sequence_file_in_f, file_view)
{
    sequence_file_in fin{std::istringstream{input}, sequence_file_format_fasta{}};

    auto minimum_length_filter = view::filter([] (auto const & rec)
    {
        return size(get<field::SEQ>(rec)) >= 5;
    });

    size_t counter = 1; // the first record will be filtered out
    for (auto & rec : fin | minimum_length_filter)
    {
        EXPECT_TRUE((ranges::equal(get<field::SEQ>(rec), seq_comp[counter])));
        EXPECT_TRUE((ranges::equal(get<field::ID>(rec),  id_comp[counter])));
        EXPECT_TRUE(empty(get<field::QUAL>(rec)));

        counter++;
    }

    EXPECT_EQ(counter, 3u);
}

TEST_F(sequence_file_in_f, column_reading)
{
    sequence_file_in fin{std::istringstream{input}, sequence_file_format_fasta{}};

    auto & seqs  = get<field::SEQ>(fin);                                    // by field
    auto & ids   = get<1>(fin);                                             // by index
    auto & quals = get<typename decltype(fin)::quality_column_type>(fin);   // by type

    ASSERT_EQ(seqs.size(), 3ul);
    ASSERT_EQ(ids.size(), 3ul);
    ASSERT_EQ(quals.size(), 3ul);

    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_TRUE((ranges::equal(seqs[i], seq_comp[i])));
        EXPECT_TRUE((ranges::equal(ids[i],  id_comp[i])));
        EXPECT_TRUE(empty(quals[i]));
    }
}

TEST_F(sequence_file_in_f, column_reading_temporary)
{
    sequence_file_in{std::istringstream{input}, sequence_file_format_fasta{}};

    auto seqs = get<field::SEQ>(sequence_file_in{std::istringstream{input}, sequence_file_format_fasta{}});

    ASSERT_EQ(seqs.size(), 3ul);

    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_TRUE((ranges::equal(seqs[i], seq_comp[i])));
    }
}

TEST_F(sequence_file_in_f, column_reading_decomposed)
{
    sequence_file_in fin{std::istringstream{input}, sequence_file_format_fasta{}};

    auto & [ seqs, ids , quals ] = fin;

    ASSERT_EQ(seqs.size(), 3ul);
    ASSERT_EQ(ids.size(), 3ul);
    ASSERT_EQ(quals.size(), 3ul);

    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_TRUE((ranges::equal(seqs[i], seq_comp[i])));
        EXPECT_TRUE((ranges::equal(ids[i],  id_comp[i])));
        EXPECT_TRUE(empty(quals[i]));
    }
}

TEST_F(sequence_file_in_f, column_reading_decomposed_temporary)
{
    auto && [ seqs, ids , quals ] = sequence_file_in{std::istringstream{input}, sequence_file_format_fasta{}};

    ASSERT_EQ(seqs.size(), 3ul);
    ASSERT_EQ(ids.size(), 3ul);
    ASSERT_EQ(quals.size(), 3ul);

    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_TRUE((ranges::equal(seqs[i], seq_comp[i])));
        EXPECT_TRUE((ranges::equal(ids[i],  id_comp[i])));
        EXPECT_TRUE(empty(quals[i]));
    }
}
