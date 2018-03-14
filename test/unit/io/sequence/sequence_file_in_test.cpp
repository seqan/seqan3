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

#include <seqan3/core/concept/iterator.hpp>
#include <seqan3/io/sequence/sequence_file_in.hpp>
#include <seqan3/range/view/to_char.hpp>

using namespace seqan3;
using namespace seqan3::literal;

TEST(sequence_file_in_iterator, concepts)
{
    using it_t = typename sequence_file_in<>::iterator;
    using sen_t = typename sequence_file_in<>::sentinel;

    EXPECT_TRUE((seqan3::input_iterator_concept<it_t>));
    EXPECT_TRUE((seqan3::sentinel_concept<sen_t, it_t>));

    // const_iterator is not a an iterator, i.e. files are not const-iterable:
    using cit_t = typename sequence_file_in<>::const_iterator;
    EXPECT_FALSE((seqan3::input_iterator_concept<cit_t>));
}


struct sequence_file_in_f : public ::testing::Test
{
    std::string input
    {
        "> TEST1\n"
        "ACGT\n"
        "> Test2\n"
        "AGGCTGN\n"
        "> Test3\n"
        "GGAGTATAATATATATATATATAT\n"
    };

    std::string seq_comp[3]
    {
        "ACGT",
        "AGGCTGN",
        "GGAGTATAATATATATATATATAT"
    };

    std::string id_comp[3]
    {
        "> TEST1",
        "> Test2",
        "> Test3"
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

TEST_F(sequence_file_in_f, record_reading_copying)
{
    /* record based reading */
    sequence_file_in fin{std::istringstream{input}, sequence_file_format_fasta{}};

    size_t counter = 0;
    for (auto [ seq, id ] : fin)
    {
        EXPECT_EQ(std::string(seq| view::to_char), seq_comp[counter]);
        EXPECT_EQ(std::string(id),                 id_comp[counter]);

        // memory still in buffer
        EXPECT_EQ(std::string(get<0>(fin.back()) | view::to_char), seq_comp[counter]);
        EXPECT_EQ(std::string(get<1>(fin.back())),                 id_comp[counter]);
        counter++;
    }
}

TEST_F(sequence_file_in_f, record_reading_ref)
{
    /* record based reading */
    sequence_file_in fin{std::istringstream{input}, sequence_file_format_fasta{}};

    size_t counter = 0;
    for (auto & [ seq, id ] : fin)
    {
        EXPECT_EQ(std::string(seq | view::to_char), seq_comp[counter]);
        EXPECT_EQ(std::string(id),                  id_comp[counter]);

        // memory still in buffer
        EXPECT_EQ(std::string(get<0>(fin.back()) | view::to_char), seq_comp[counter]);
        EXPECT_EQ(std::string(get<1>(fin.back())),                 id_comp[counter]);
        counter++;
    }
}

TEST_F(sequence_file_in_f, record_reading_move)
{
    /* record based reading */
    sequence_file_in fin{std::istringstream{input}, sequence_file_format_fasta{}};

    size_t counter = 0;
    for (auto [ seq, id ] : std::move(fin))
    {
        EXPECT_EQ(std::string(seq | view::to_char), seq_comp[counter]);
        EXPECT_EQ(std::string(id),                  id_comp[counter]);
        counter++;

        // buffers cleared TODO doesnt work
//         EXPECT_TRUE(empty(std::string(get<0>(fin.back()) | view::to_char)));
//         EXPECT_TRUE(empty(std::string(get<1>(fin.back()))));
    }
}

TEST_F(sequence_file_in_f, column_reading)
{
    /* column based reading */
    sequence_file_in fin2{std::istringstream{input}, sequence_file_format_fasta{}};

    auto [ seqs, ids ] = std::move(fin2);

    EXPECT_EQ(seqs.size(), 3ul);
    EXPECT_EQ(ids.size(), 3ul);

    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_EQ(std::string(seqs[i] | view::to_char), seq_comp[i]);
        EXPECT_EQ(std::string(ids[i]),                  id_comp[i]);
    }

/*
    // can't read a second time
    auto [ seqs2, ids2 ] = std::move(fin2);

    EXPECT_TRUE(empty(seqs2));
    EXPECT_TRUE(empty(ids2));*/
}
