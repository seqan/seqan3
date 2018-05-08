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

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence/sequence_file_format_concept.hpp>
#include <seqan3/io/sequence/sequence_file_format_fasta.hpp>
#include <seqan3/range/view/convert.hpp>

using namespace seqan3;
using namespace seqan3::literal;

TEST(general, concepts)
{
    EXPECT_TRUE((sequence_file_format_concept<sequence_file_format_fasta>));
}

struct read : public ::testing::Test
{
    std::vector<std::string> expected_ids
    {
        { "ID1" },
        { "ID2" },
        { "ID3 lala" },
    };

    std::vector<dna5_vector> expected_seqs
    {
        { "ACGTTTTTTTTTTTTTTT"_dna5 },
        { "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"_dna5 },
        { "ACGTTTA"_dna5 },
    };

    sequence_file_format_fasta format;

    detail::sequence_file_in_options_type_dummy options;

    std::string id;
    dna5_vector seq;

    void do_read_test(std::string const & input)
    {
        std::stringstream istream{input};

        for (unsigned i = 0; i < 3; ++i)
        {
            id.clear();
            seq.clear();

            EXPECT_NO_THROW(( format.read(istream, options, seq, id, std::ignore, std::ignore) ));

            EXPECT_TRUE((ranges::equal(seq, expected_seqs[i])));
            EXPECT_TRUE((ranges::equal(id, expected_ids[i])));
        }
    }
};

TEST_F(read, standard)
{
    std::string input
    {
        "> ID1\n"
        "ACGTTTTTTTTTTTTTTT\n"
        "> ID2\n"
        "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        "> ID3 lala\n"
        "ACGTTTA\n"
    };

    do_read_test(input);
}

TEST_F(read, newline_before_eof)
{
    std::string input
    {
        "> ID1\n"
        "ACGTTTTTTTTTTTTTTT\n"
        "> ID2\n"
        "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        "> ID3 lala\n"
        "ACGTTTA"
    };
    do_read_test(input);
}

TEST_F(read, noblank_before_id)
{
    std::string input
    {
        ">ID1\n"
        "ACGTTTTTTTTTTTTTTT\n"
        ">ID2\n"
        "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        ">ID3 lala\n"
        "ACGTTTA\n"
    };
    do_read_test(input);
}

TEST_F(read, whitespace_in_seq)
{
    std::string input
    {
        "> ID1\n"
        "ACGTTTT\n\nTTTTTTTTTTT\n"
        "\n"
        "> ID2\n"
        "ACGTTTT\t\tTTTTTTTTTTT\t\nTTTTTTTTTTT\vTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\rTTTTTTTTTTTTTTTTT\n"
        "> ID3 lala\n"
        "ACGT\fTTA\n"
    };
    do_read_test(input);
}

TEST_F(read, digits_in_seq)
{
    std::string input
    {
        "> ID1\n"
        "10  ACGTTTTTTTTTTTTTTT\n"
        "> ID2\n"
        "  80 ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT  900"
        "1000 TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        "> ID3 lala\n"
        "ACGT9T5T2A\n"
    };
    do_read_test(input);
}

TEST_F(read, old_id_style)
{
    std::string input
    {
        "; ID1\n"
        "ACGTTTTTTTTTTTTTTT\n"
        "; ID2\n"
        "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        "; ID3 lala\n"
        "ACGTTTA\n"
    };

    do_read_test(input);
}

TEST_F(read, mixed_issues)
{
    std::string input
    {
        "> ID1\n"
        "ACGTTTT\n\nTTTTTTTTTTT\n"
        "\n"
        ";ID2\n"
        "ACGTTTT\t75\tTTTTTTTTTTT\t\nTTTTTTTTTTT9\vTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\rTTTTTTTTTTTTTTTTT\n"
        "> ID3 lala\n"
        "ACGT\fTTA"
    };
    do_read_test(input);
}

TEST_F(read, options_truncate_ids)
{
    std::string input
    {
        "> ID1\n"
        "ACGTTTTTTTTTTTTTTT\n"
        "> ID2\n"
        "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        "> ID3 lala\n"
        "ACGTTTA\n"
    };

    options.truncate_ids = true;
    expected_ids[2] = "ID3"; // "lala" is stripped
    do_read_test(input);
}

TEST_F(read, only_seq)
{
    std::string input
    {
        "> ID1\n"
        "ACGTTTTTTTTTTTTTTT\n"
        "> ID2\n"
        "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        "> ID3 lala\n"
        "ACGTTTA\n"
    };

    std::stringstream istream{input};

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq.clear();

        format.read(istream, options, seq, std::ignore, std::ignore, std::ignore);

        EXPECT_TRUE((ranges::equal(seq, expected_seqs[i])));
    }
}

TEST_F(read, only_id)
{
    std::string input
    {
        "> ID1\n"
        "ACGTTTTTTTTTTTTTTT\n"
        "> ID2\n"
        "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        "> ID3 lala\n"
        "ACGTTTA\n"
    };

    std::stringstream istream{input};

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq.clear();

        format.read(istream, options, std::ignore, id, std::ignore, std::ignore);

        EXPECT_TRUE((ranges::equal(id, expected_ids[i])));
    }
}

TEST_F(read, seq_qual)
{
    std::string input
    {
        "> ID1\n"
        "ACGTTTTTTTTTTTTTTT\n"
        "> ID2\n"
        "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        "> ID3 lala\n"
        "ACGTTTA\n"
    };

    std::stringstream istream{input};

    std::vector<quality_composition<dna5, illumina18>> seq_qual;

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq_qual.clear();

        format.read(istream, options, std::ignore, id, std::ignore, seq_qual);

        EXPECT_TRUE((ranges::equal(id, expected_ids[i])));
        EXPECT_TRUE((ranges::equal(seq_qual | view::convert<dna5>, expected_seqs[i])));
    }
}

TEST_F(read, fail_no_id)
{
    std::string input
    {
        "! ID1\n"
        "ACGTTTTTTTTTTTTTTT\n"
        "> ID2\n"
        "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        "> ID3 lala\n"
        "ACGTTTA\n"
    };

    std::stringstream istream{input};

    EXPECT_THROW( (format.read(istream, options, std::ignore, std::ignore, std::ignore, std::ignore)),
                  parse_error );
}

TEST_F(read, fail_wrong_char)
{
    std::string input
    {
        "> ID1\n"
        "ACGPTTTTTTTTTTTTTT\n"
        "> ID2\n"
        "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        "> ID3 lala\n"
        "ACGTTTA\n"
    };

    std::stringstream istream{input};

    EXPECT_THROW( (format.read(istream, options, seq, id, std::ignore, std::ignore)),
                  parse_error );
}
