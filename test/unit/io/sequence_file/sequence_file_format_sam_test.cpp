// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/test/pretty_printing.hpp>

#include "sequence_file_format_test_template.hpp"

template <>
struct sequence_file_read<seqan3::format_sam> : public sequence_file_data
{
    std::string standard_input
    {
R"(ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*+,-./++-
ID2	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE
ID3 lala	0	*	0	0	*	*	0	0	ACGTTTA	!!!!!!!
)"
    };

    std::string illegal_alphabet_character_input
    {
R"(ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTT?TTTTTTT	!##$%&'()*+,-./++-
)"
    };

    std::string standard_output
    {
R"(ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*+,-./++-
ID2	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE
ID3 lala	0	*	0	0	*	*	0	0	ACGTTTA	!!!!!!!
)"
    };

    std::string no_or_ill_formatted_id_input
    {
R"(*	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*+,-./++-
)"
	};
};

// ---------------------------------------------------------------------------------------------------------------------
// parametrized tests
// ---------------------------------------------------------------------------------------------------------------------

INSTANTIATE_TYPED_TEST_SUITE_P(sam, sequence_file_read, seqan3::format_sam, );
INSTANTIATE_TYPED_TEST_SUITE_P(sam, sequence_file_write, seqan3::format_sam, );

// ----------------------------------------------------------------------------
// reading
// ----------------------------------------------------------------------------
struct read_sam : sequence_file_read<seqan3::format_sam>
{
    seqan3::sequence_file_input_options<seqan3::dna5, false> options{};

    void do_read_test(std::string const & input)
    {
        std::stringstream istream{input};
        seqan3::sequence_file_input fin{istream, seqan3::format_sam{}};

        auto it = fin.begin();
        for (unsigned i = 0; i < 3; ++i, it++)
        {
            EXPECT_EQ(seqan3::get<seqan3::field::seq>(*it), seqs[i]);
            EXPECT_EQ(seqan3::get<seqan3::field::id>(*it), ids[i]);
            EXPECT_EQ(seqan3::get<seqan3::field::qual>(*it), quals[i]);
        }
    }
};

TEST_F(read_sam, tags)
{
    std::string input
    {
R"(ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*+,-./++-	FI:i:1
ID2	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE	AS:i:3
ID3 lala	0	*	0	0	*	*	0	0	ACGTTTA	!!!!!!!	TI:i:2
)"
    };
    do_read_test(input);
}

TEST_F(read_sam, mixed_issues)
{
  std::string input
  {
R"(ID1	0	BABABA	200	0	*	BABABA	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*+,-./++-	FI:i:1
ID2	0	*	0	0	BABA	*	30	0	ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE
ID3 lala	0	*	0	0	HAHAHAHA+	*	0	0	ACGTTTA	!!!!!!!
)"
    };
    do_read_test(input);
}

TEST_F(read_sam, no_qual)
{
   std::string input
   {
R"(ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	*
ID2	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDE
ID3 lala	0	*	0	0	*	*	0	0	ACGTTTA	!!!!!!!
)"
    };

    quals[0] = { ""_phred42 };
    do_read_test(input);
}

TEST_F(read_sam, qual_too_short)
{
    std::stringstream istream{std::string{
R"(ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()-./++-
)"}};

    seqan3::sequence_file_input fin{istream, seqan3::format_sam{}};
    EXPECT_THROW(fin.begin(), seqan3::format_error);
}

TEST_F(read_sam, qual_too_long)
{
    std::stringstream istream{std::string{
R"(ID1	0	*	0	0	*	*	0	0	ACGTTTTTTTTTTTTTTT	!##$%&'()*+,-./++-+
)"}};

    seqan3::sequence_file_input fin{istream, seqan3::format_sam{}};
    EXPECT_THROW(fin.begin(), seqan3::format_error);
}

TEST_F(read_sam, no_seq)
{
    std::stringstream istream{std::string{
R"(ID 1	0	*	0	0	*	*	0	0	*	!##$%&'()*+,-./++-
)"}};

    seqan3::sequence_file_input fin{istream, seqan3::format_sam{}};
    EXPECT_THROW(fin.begin(), seqan3::parse_error);
}

// ----------------------------------------------------------------------------
// writing
// ----------------------------------------------------------------------------
struct write : public ::testing::Test
{
    std::vector<seqan3::dna5_vector> seqs
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

    std::vector<std::vector<seqan3::phred42>> quals
    {
        { "!##$"_phred42 },
        { "!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE"_phred42 },
        { "!!*+,-./+*+,-./+!!FF!!"_phred42 },
    };

    seqan3::sequence_file_output_options options{};
    std::ostringstream ostream;

    void do_write_test()
    {
    	seqan3::sequence_file_output fout{ostream, seqan3::format_sam{}};
        for (unsigned i = 0; i < 3; ++i)
            EXPECT_NO_THROW(( fout.emplace_back(seqs[i], ids[i], quals[i]) ));
        ostream.flush();
    }

    void do_write_test_no_qual()
    {
        seqan3::sequence_file_output fout{ostream, seqan3::format_sam{}};
        for (unsigned i = 0; i < 3; ++i)
            EXPECT_NO_THROW(( fout.emplace_back(seqs[i], ids[i], ""_phred42) ));
        ostream.flush();
    }
};

TEST_F(write, no_id)
{
    std::string comp
    {
R"(*	0	*	0	0	*	*	0	0	ACGT	!##$
Test2	0	*	0	0	*	*	0	0	AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN	!##$&'()*+,-./+)*+,-)*+,-)*+,-)*+,BDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDEBDBDDEBDBEEBEBE
Test3	0	*	0	0	*	*	0	0	GGAGTATAATATATATATATAT	!!*+,-./+*+,-./+!!FF!!
)"
    };

    ids[0] = "";
    do_write_test();

    EXPECT_EQ(ostream.str(), comp);
}

TEST_F(write, no_seq)
{
    std::string comp
    {
R"(TEST 1	0	*	0	0	*	*	0	0	*	*
Test2	0	*	0	0	*	*	0	0	AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN	*
Test3	0	*	0	0	*	*	0	0	GGAGTATAATATATATATATAT	*
)"
    };

    seqs[0] = ""_dna5;
    do_write_test_no_qual();

    EXPECT_EQ(ostream.str(), comp);
}

// No qualities given
TEST_F(write, no_qual)
{
    std::string comp
    {
R"(TEST 1	0	*	0	0	*	*	0	0	ACGT	*
Test2	0	*	0	0	*	*	0	0	AGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN	*
Test3	0	*	0	0	*	*	0	0	GGAGTATAATATATATATATAT	*
)"
    };

    do_write_test_no_qual();

    EXPECT_EQ(ostream.str(), comp);
}
