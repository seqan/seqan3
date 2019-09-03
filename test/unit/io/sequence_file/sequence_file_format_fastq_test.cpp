// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/input_format_concept.hpp>
#include <seqan3/io/sequence_file/output_format_concept.hpp>
#include <seqan3/io/sequence_file/format_fastq.hpp>
#include <seqan3/range/view/convert.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

using namespace seqan3;

// ----------------------------------------------------------------------------
// general
// ----------------------------------------------------------------------------

TEST(general, concepts)
{
    EXPECT_TRUE((sequence_file_input_format<format_fastq>));
    EXPECT_TRUE((sequence_file_output_format<format_fastq>));
}

// ----------------------------------------------------------------------------
// reading
// ----------------------------------------------------------------------------

struct read : public ::testing::Test
{
    std::string input
    {
        "@ID1\n"
        "ACGTTTTTTTTTTTTTTT\n"
        "+\n"
        "!##$%&'()*+,-./++-\n"
        "@ID2\n"
        "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        "+\n"
        "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
        "@ID3 lala\n"
        "ACGTTTA\n"
        "+\n"
        "!!!!!!!\n"
    };

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

    std::vector<std::vector<phred42>> expected_quals
    {
        { "!##$%&'()*+,-./++-"_phred42 },
        { "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE"_phred42 },
        { "!!!!!!!"_phred42 },
    };

    detail::sequence_file_input_format_REMOVEME<format_fastq> format;

    sequence_file_input_options<dna5, false> options;

    std::string id;
    dna5_vector seq;
    std::vector<phred42> qual;

    void do_read_test(std::string const & input)
    {
        std::stringstream istream{input};

        for (unsigned i = 0; i < 3; ++i)
        {
            id.clear();
            seq.clear();
            qual.clear();

            EXPECT_NO_THROW(( format.read(istream, options, seq, id, qual) ));

            EXPECT_TRUE((std::ranges::equal(seq, expected_seqs[i])));
            EXPECT_TRUE((std::ranges::equal(id, expected_ids[i])));
            EXPECT_TRUE((std::ranges::equal(qual, expected_quals[i])));
        }
    }
};

TEST_F(read, standard)
{
    do_read_test(input);
}

TEST_F(read, newline_before_eof)
{
    input =
    {
        "@ID1\n"
        "ACGTTTTTTTTTTTTTTT\n"
        "+\n"
        "!##$%&'()*+,-./++-\n"
        "@ID2\n"
        "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        "+\n"
        "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
        "@ID3 lala\n"
        "ACGTTTA\n"
        "+\n"
        "!!!!!!!"
    };
    do_read_test(input);
}

TEST_F(read, whitespace_in_seq_qual)
{
    input =
    {
        "@ID1\n"
        "ACGTTTTTTTT\nTTTTTTT\n"
        "+\n"
        "!##$\n%&'()*+,-./++-\n"
        "@ID2\n"
        "ACGTTTTTTTTTT\r\nTTTTTT\r\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\r\nTTTTTTTTTTTTTTTTTTTTTTT\n"
        "+\n"
        "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBD\nEBDEB\nDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
        "@ID3 lala\n"
        "ACGTT\nTA\n"
        "+\n"
        "!!!!!\n!!\n"
    };
    do_read_test(input);
}

TEST_F(read, double_id_style)
{
    input =
    {
        "@ID1\n"
        "ACGTTTTTTTTTTTTTTT\n"
        "+ID1\n"
        "!##$%&'()*+,-./++-\n"
        "@ID2\n"
        "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        "+ID2\n"
        "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
        "@ID3 lala\n"
        "ACGTTTA\n"
        "+ID3 lala\n"
        "!!!!!!!\n"
    };
    do_read_test(input);
}

TEST_F(read, mixed_issues)
{
    input =
    {
        "@ID1\n"
        "ACGTTTTTTTT\nTTTTTTT\n"
        "+\n"
        "!##$\n%&'()*+,-./++-\n"
        "@ID2\n"
        "ACGTTTTTTTTTT\r\nTTTTTT\r\nTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\r\nTTTTTTTTTTTTTTTTTTTTTTT\n"
        "+\n"
        "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBD\nEBDEB\nDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
        "@ID3 lala\n"
        "ACGTT\nTA\n"
        "+ID3 lala\n"
        "!!!!!\n!!"
    };
    do_read_test(input);
}

TEST_F(read, options_truncate_ids)
{
    input =
    {
        "@ID1\n"
        "ACGTTTTTTTTTTTTTTT\n"
        "+\n"
        "!##$%&'()*+,-./++-\n"
        "@ID2\n"
        "ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n"
        "+\n"
        "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE\n"
        "@ID3 lala\n"
        "ACGTTTA\n"
        "+\n"
        "!!!!!!!\n"
    };

    options.truncate_ids = true;
    expected_ids[2] = "ID3"; // "lala" is stripped
    do_read_test(input);
}

TEST_F(read, only_seq)
{
    std::stringstream istream{input};

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq.clear();
        qual.clear();

        format.read(istream, options, seq, std::ignore, std::ignore);

        EXPECT_TRUE((std::ranges::equal(seq, expected_seqs[i])));
    }
}

TEST_F(read, only_id)
{
    std::stringstream istream{input};

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq.clear();
        qual.clear();

        format.read(istream, options, std::ignore, id, std::ignore);

        EXPECT_TRUE((std::ranges::equal(id, expected_ids[i])));
    }
}

TEST_F(read, only_qual)
{
    std::stringstream istream{input};

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq.clear();
        qual.clear();

        format.read(istream, options, std::ignore, std::ignore, qual);

        EXPECT_TRUE((std::ranges::equal(qual, expected_quals[i])));
    }
}

TEST_F(read, seq_qual)
{
    std::stringstream istream{input};

    std::vector<qualified<dna5, phred42>> seq_qual;
    sequence_file_input_options<dna5, true> options2;

    for (unsigned i = 0; i < 3; ++i)
    {
        id.clear();
        seq_qual.clear();

        format.read(istream, options2, seq_qual, id, seq_qual);

        EXPECT_TRUE((std::ranges::equal(id, expected_ids[i])));
        EXPECT_TRUE((std::ranges::equal(seq_qual | view::convert<dna5>, expected_seqs[i])));
        EXPECT_TRUE((std::ranges::equal(seq_qual | view::convert<phred42>, expected_quals[i])));
    }
}

TEST_F(read, fail_no_id)
{
    input[0] = '#'; // invalid character for ID line

    std::stringstream istream{input};

    EXPECT_THROW( (format.read(istream, options, std::ignore, id, std::ignore)),
                  parse_error );
}

//TODO fail_no_2nd_id
TEST_F(read, fail_no_seq_after_id)
{
    input =
    {
        "@ID1\n"
        "ACGTTTTTTTTTTTTTTT"
    };

    std::stringstream istream{input};

    EXPECT_THROW( (format.read(istream, options, seq, id, qual)),
                  unexpected_end_of_input );
}

//TODO fail_no_quals
//TODO fail_quals_shorter_seq

TEST_F(read, fail_wrong_char)
{
    input =
    {
        "@ID1\n"
        "ACGTTPTTTTTTTTTTTT\n"
        "+\n"
        "!##$%&'()*+,-./++-\n"
    };

    std::stringstream istream{input};

    EXPECT_THROW( (format.read(istream, options, seq, id, qual)),
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
        "GGAGTATAATATATATATATAT"_dna5
    };

    std::vector<std::string> ids
    {
        "TEST 1",
        "Test2",
        "Test3"
    };

    std::vector<std::vector<phred42>> quals
    {
        { "!##$"_phred42 },
        { "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE"_phred42 },
        { "!!*+,-./+*+,-./+!!FF!!"_phred42 },
    };

    detail::sequence_file_output_format_REMOVEME<format_fastq> format;

    sequence_file_output_options options;

    std::ostringstream ostream;

    void do_write_test()
    {
        for (unsigned i = 0; i < 3; ++i)
            EXPECT_NO_THROW(( format.write(ostream, options, seqs[i], ids[i], quals[i]) ));

        ostream.flush();
    }
};

TEST_F(write, arg_handling_id_missing)
{
    EXPECT_THROW( (format.write(ostream, options, seqs[0], std::ignore, quals[0])),
                   std::logic_error );
}

TEST_F(write, arg_handling_id_empty)
{
    EXPECT_THROW( (format.write(ostream, options, seqs[0], std::string_view{""}, quals[0])),
                   std::runtime_error );
}

TEST_F(write, arg_handling_seq_missing)
{
    EXPECT_THROW( (format.write(ostream, options, std::ignore, ids[0], quals[0])),
                   std::logic_error );
}

TEST_F(write, arg_handling_seq_empty)
{
    EXPECT_THROW( (format.write(ostream, options, std::string_view{""}, ids[0], quals[0])),
                   std::runtime_error );
}

TEST_F(write, arg_handling_qual_missing)
{
    EXPECT_THROW( (format.write(ostream, options, seqs[0], ids[0], std::ignore)),
                   std::logic_error );
}

TEST_F(write, arg_handling_qual_empty)
{
    EXPECT_THROW( (format.write(ostream, options, seqs[0], ids[0], std::string_view{""})),
                   std::runtime_error );
}

TEST_F(write, default_options)
{
    std::string comp
    {
        "@TEST 1\n"
        "ACGT\n"
        "+\n"
        "!##$\n"
        "@Test2\n"
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN\n"
        "+\n"
        "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE\n"
        "@Test3\n"
        "GGAGTATAATATATATATATAT\n"
        "+\n"
        "!!*+,-./+*+,-./+!!FF!!\n"
    };

    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, seq_qual)
{
    std::vector<std::vector<qualified<dna5, phred42>>> vec;
    vec.resize(3);
    for (unsigned i = 0; i < 3; ++i)
    {
        vec[i].resize(seqs[i].size());
        for (unsigned j = 0; j < seqs[i].size(); ++j)
        {
            vec[i][j] = seqs[i][j];
            vec[i][j] = quals[i][j];
        }
    }

    for (unsigned i = 0; i < 3; ++i)
        EXPECT_NO_THROW(( format.write(ostream,
                                       options,
                                       vec[i] | view::convert<dna5>,
                                       ids[i],
                                       vec[i] | view::convert<phred42>) ));

    ostream.flush();

    std::string comp
    {
        "@TEST 1\n"
        "ACGT\n"
        "+\n"
        "!##$\n"
        "@Test2\n"
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN\n"
        "+\n"
        "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE\n"
        "@Test3\n"
        "GGAGTATAATATATATATATAT\n"
        "+\n"
        "!!*+,-./+*+,-./+!!FF!!\n"
    };

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, options_fastq_double_id)
{
    options.fastq_double_id = true;

    std::string comp
    {
        "@TEST 1\n"
        "ACGT\n"
        "+TEST 1\n"
        "!##$\n"
        "@Test2\n"
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN\n"
        "+Test2\n"
        "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE\n"
        "@Test3\n"
        "GGAGTATAATATATATATATAT\n"
        "+Test3\n"
        "!!*+,-./+*+,-./+!!FF!!\n"
    };

    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, options_add_carriage_return)
{
    options.add_carriage_return = true;

    std::string comp
    {
        "@TEST 1\r\n"
        "ACGT\r\n"
        "+\r\n"
        "!##$\r\n"
        "@Test2\r\n"
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN\r\n"
        "+\r\n"
        "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE\r\n"
        "@Test3\r\n"
        "GGAGTATAATATATATATATAT\r\n"
        "+\r\n"
        "!!*+,-./+*+,-./+!!FF!!\r\n"
    };

    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, options_all)
{
    options.add_carriage_return = true;
    options.fastq_double_id = true;

    std::string comp
    {
        "@TEST 1\r\n"
        "ACGT\r\n"
        "+TEST 1\r\n"
        "!##$\r\n"
        "@Test2\r\n"
        "AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN\r\n"
        "+Test2\r\n"
        "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE\r\n"
        "@Test3\r\n"
        "GGAGTATAATATATATATATAT\r\n"
        "+Test3\r\n"
        "!!*+,-./+*+,-./+!!FF!!\r\n"
    };
    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}
