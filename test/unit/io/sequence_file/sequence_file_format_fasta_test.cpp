// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>
#include <range/v3/view/transform.hpp>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/io/sequence_file/format_fasta.hpp>
#include <seqan3/range/view/convert.hpp>

using namespace seqan3;

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TEST(general, concepts)
{
    EXPECT_TRUE((sequence_file_input_format_concept<sequence_file_format_fasta>));
    EXPECT_TRUE((sequence_file_output_format_concept<sequence_file_format_fasta>));
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

    sequence_file_input_options<dna5, false> options;

    std::string id;
    dna5_vector seq;

    void do_read_test(std::string const & input)
    {
        std::stringstream istream{input};

        for (unsigned i = 0; i < 3; ++i)
        {
            id.clear();
            seq.clear();

            EXPECT_NO_THROW(( format.read(istream, options, seq, id, std::ignore) ));

            EXPECT_TRUE((std::ranges::equal(seq, expected_seqs[i])));
            EXPECT_TRUE((std::ranges::equal(id, expected_ids[i])));
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

        format.read(istream, options, seq, std::ignore, std::ignore);

        EXPECT_TRUE((std::ranges::equal(seq, expected_seqs[i])));
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

        format.read(istream, options, std::ignore, id, std::ignore);

        EXPECT_TRUE((std::ranges::equal(id, expected_ids[i])));
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
    sequence_file_input_options<dna5, true> options2;

    std::vector<qualified<dna5, phred42>> seq_qual;

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq_qual.clear();

        format.read(istream, options2, seq_qual, id, seq_qual);

        EXPECT_TRUE((std::ranges::equal(id, expected_ids[i])));
        EXPECT_TRUE((std::ranges::equal(seq_qual | view::convert<dna5>, expected_seqs[i])));
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

    EXPECT_THROW( (format.read(istream, options, std::ignore, std::ignore, std::ignore)),
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

    EXPECT_THROW( (format.read(istream, options, seq, id, std::ignore)),
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

    sequence_file_output_options options;

    std::ostringstream ostream;

    void do_write_test()
    {
        for (unsigned i = 0; i < 3; ++i)
            EXPECT_NO_THROW(( format.write(ostream, options, seqs[i], ids[i], std::ignore) ));

        ostream.flush();
    }
};

TEST_F(write, arg_handling_id_missing)
{
    EXPECT_THROW( (format.write(ostream, options, seqs[0], std::ignore, std::ignore)),
                   std::logic_error );
}

TEST_F(write, arg_handling_id_empty)
{
    EXPECT_THROW( (format.write(ostream, options, seqs[0], std::string_view{""}, std::ignore)),
                   std::runtime_error );
}

TEST_F(write, arg_handling_seq_missing)
{
    EXPECT_THROW( (format.write(ostream, options, std::ignore, ids[0], std::ignore)),
                   std::logic_error );
}

TEST_F(write, arg_handling_seq_empty)
{
    EXPECT_THROW( (format.write(ostream, options, std::string_view{""}, ids[0], std::ignore)),
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

    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, seq_qual)
{
    auto convert_to_qualified = ranges::view::transform([] (auto const in)
    {
        return qualified<dna5, phred42>{} = in;
    });

    for (unsigned i = 0; i < 3; ++i)
        EXPECT_NO_THROW(( format.write(ostream,
                                       options,
                                       seqs[i] | convert_to_qualified,
                                       ids[i],
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

    do_write_test();

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

    do_write_test();

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

    do_write_test();

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

    do_write_test();

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
    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}
