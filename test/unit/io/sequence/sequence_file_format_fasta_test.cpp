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
#include <range/v3/view/transform.hpp>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence/sequence_file_in_format_concept.hpp>
#include <seqan3/io/sequence/sequence_file_out_format_concept.hpp>
#include <seqan3/io/sequence/sequence_file_format_fasta.hpp>
#include <seqan3/range/view/convert.hpp>

using namespace seqan3;
using namespace seqan3::literal;

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TEST(general, concepts)
{
    EXPECT_TRUE((sequence_file_in_format_concept<sequence_file_format_fasta>));
    EXPECT_TRUE((sequence_file_out_format_concept<sequence_file_format_fasta>));
}

// ----------------------------------------------------------------------------
// reading
// ----------------------------------------------------------------------------

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

    sequence_file_in_options<dna5> options;

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

// ----------------------------------------------------------------------------
// writing
// ----------------------------------------------------------------------------

struct write : public ::testing::Test
{
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

    sequence_file_format_fasta format;

    sequence_file_out_options options;

    std::ostringstream ostream;

    void do_read_test()
    {
        for (unsigned i = 0; i < 3; ++i)
            EXPECT_NO_THROW(( format.write(ostream, options, seqs[i], ids[i], std::ignore, std::ignore) ));

        ostream.flush();
    }
};

TEST_F(write, arg_handling_id_missing)
{
    EXPECT_THROW( (format.write(ostream, options, seqs[0], std::ignore, std::ignore, std::ignore)),
                   std::logic_error );
}

TEST_F(write, arg_handling_id_empty)
{
    EXPECT_THROW( (format.write(ostream, options, seqs[0], std::string_view{""}, std::ignore, std::ignore)),
                   std::runtime_error );
}

TEST_F(write, arg_handling_seq_missing)
{
    EXPECT_THROW( (format.write(ostream, options, std::ignore, ids[0], std::ignore, std::ignore)),
                   std::logic_error );
}

TEST_F(write, arg_handling_seq_empty)
{
    EXPECT_THROW( (format.write(ostream, options, std::string_view{""}, ids[0], std::ignore, std::ignore)),
                   std::runtime_error );
}

TEST_F(write, arg_handling_seq_qual_empty)
{
    EXPECT_THROW( (format.write(ostream, options, std::ignore, ids[0], std::ignore, std::string_view{""})),
                   std::runtime_error );
}

TEST_F(write, default_options)
{
    std::string comp
    {
        "> TEST 1\n"
        "ACGT\n"
        "> Test2\n"
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGG\nCTGNAGGCTGN\n"
        //                                             linebreak inserted after 80 char  ^
        "> Test3\n"
        "GGAGTATAATATATATATATATAT\n"
    };

    do_read_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, seq_qual)
{
    auto convert_to_qualified = ranges::view::transform([] (auto const in)
    {
        return quality_composition<dna5, illumina18>{} = in;
    });

    for (unsigned i = 0; i < 3; ++i)
        EXPECT_NO_THROW(( format.write(ostream,
                                    options,
                                    std::ignore,
                                    ids[i],
                                    std::ignore,
                                    seqs[i] | convert_to_qualified) ));

    ostream.flush();

    std::string comp
    {
        "> TEST 1\n"
        "ACGT\n"
        "> Test2\n"
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGG\nCTGNAGGCTGN\n"
        //                                             linebreak inserted after 80 char  ^
        "> Test3\n"
        "GGAGTATAATATATATATATATAT\n"
    };

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, options_letters_per_line)
{
    options.fasta_letters_per_line = 7;

    std::string comp
    {
        "> TEST 1\n"
        "ACGT\n"
        "> Test2\n"
        "AGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\nAGGCTGN\n"
        "AGGCTGN\n"
        "> Test3\n"
        "GGAGTAT\nAATATAT\nATATATA\nTAT\n"
    };

    do_read_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, options_legacy_id_marker)
{
    options.fasta_legacy_id_marker = true;

    std::string comp
    {
        "; TEST 1\n"
        "ACGT\n"
        "; Test2\n"
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGG\nCTGNAGGCTGN\n"
        //                                             linebreak inserted after 80 char  ^
        "; Test3\n"
        "GGAGTATAATATATATATATATAT\n"
    };

    do_read_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, options_blank_before_id)
{
    options.fasta_blank_before_id = false;

    std::string comp
    {
        ">TEST 1\n"
        "ACGT\n"
        ">Test2\n"
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGG\nCTGNAGGCTGN\n"
        //                                             linebreak inserted after 80 char  ^
        ">Test3\n"
        "GGAGTATAATATATATATATATAT\n"
    };

    do_read_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, options_add_carriage_return)
{
    options.add_carriage_return = true;

    std::string comp
    {
        "> TEST 1\r\n"
        "ACGT\r\n"
        "> Test2\r\n"
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGG\r\nCTGNAGGCTGN\r\n"
        //                                             linebreak inserted after 80 char  ^
        "> Test3\r\n"
        "GGAGTATAATATATATATATATAT\r\n"
    };

    do_read_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, options_all)
{
    options.add_carriage_return = true;
    options.fasta_blank_before_id = false;
    options.fasta_legacy_id_marker = true;
    options.fasta_letters_per_line = 21;

    std::string comp
    {
        ";TEST 1\r\n"
        "ACGT\r\n"
        ";Test2\r\n"
        "AGGCTGNAGGCTGNAGGCTGN\r\nAGGCTGNAGGCTGNAGGCTGN\r\nAGGCTGNAGGCTGNAGGCTGN\r\nAGGCTGNAGGCTGNAGGCTGN\r\n"
        "AGGCTGN\r\n"
        ";Test3\r\n"
        "GGAGTATAATATATATATATA\r\nTAT\r\n"
    };
    do_read_test();

    EXPECT_EQ(ostream.str(), comp);
}
