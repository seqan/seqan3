// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/test/tmp_filename.hpp>
#include <seqan3/std/iterator>

TEST(rows, assign_sequence_files)
{
    std::string const input
    {
        ">TEST 1\n"
        "ACGT\n"
        "> Test2\n"
        "AGGCTGN AGGCTGN AGGCTGN AGGCTGN AGGCTGN AGGCTGN AGGCTGN AGGCTGN AGGCTGN AGGCTGN AGGCTGN AGGCTGN AGGCTGN\n\n"
        "> Test3\n"
        "GGAGTATAATATATATATATATAT\n"
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

    seqan3::sequence_file_input fin{std::istringstream{input}, seqan3::format_fasta{}};
    seqan3::sequence_file_output fout{std::ostringstream{}, seqan3::format_fasta{}};
    fout.options.fasta_letters_per_line = 0;

    fout = fin;

    fout.get_stream().flush();
    EXPECT_EQ(reinterpret_cast<std::ostringstream&>(fout.get_stream()).str(), output_comp);
}

TEST(integration, assign_sequence_file_pipes)
{
    std::string const input
    {
        "> TEST1\n"
        "ACGT\n"
        "> Test2\n"
        "AGGCTGNAGGCTGAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN\n"
        "> Test3\n"
        "GGAGTATAATATATATATATATAT\n"
    };

    // valid without assignment?
    seqan3::sequence_file_input{std::istringstream{input}, seqan3::format_fasta{}} |
        seqan3::sequence_file_output{std::ostringstream{}, seqan3::format_fasta{}};

    // valid with assignment and check contents
    auto fout = seqan3::sequence_file_input{std::istringstream{input}, seqan3::format_fasta{}} |
                seqan3::sequence_file_output{std::ostringstream{}, seqan3::format_fasta{}};

    fout.get_stream().flush();
    EXPECT_EQ(reinterpret_cast<std::ostringstream&>(fout.get_stream()).str(), input);
}

TEST(integration, view)
{
    std::string const input
    {
        "> TEST1\n"
        "ACGT\n"
        "> Test2\n"
        "AGGCTGNAGGCTGAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN\n"
        "> Test3\n"
        "GGAGTATAATATATATATATATAT\n"
    };

    std::string const output
    {
        "> TEST1\n"
        "ACGT\n"
        "> Test2\n"
        "AGGCTGNAGGCTGAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGNAGGCTGN\n"
    };

    // valid without assignment?
    seqan3::sequence_file_input{std::istringstream{input},
                                seqan3::format_fasta{}} | seqan3::views::persist
                                                        | seqan3::views::take(2)
                                                        | seqan3::sequence_file_output{std::ostringstream{},
                                                                                       seqan3::format_fasta{}};

    // valid with assignment and check contents
    auto fout = seqan3::sequence_file_input{std::istringstream{input}, seqan3::format_fasta{}}
              | seqan3::views::persist
              | seqan3::views::take(2)
              | seqan3::sequence_file_output{std::ostringstream{}, seqan3::format_fasta{}};

    fout.get_stream().flush();
    EXPECT_EQ(reinterpret_cast<std::ostringstream&>(fout.get_stream()).str(), output);
}

TEST(integration, convert_fastq_to_fasta)
{
    std::string const fastq_in
    {
    "@ID1\n"
    "ACGTT\n"
    "+\n"
    "!##$%\n"
    "@ID2\n"
    "TATTA\n"
    "+\n"
    ",BDEB\n"
    };

    std::string const fasta_out
    {
        "> ID1\n"
        "ACGTT\n"
        "> ID2\n"
        "TATTA\n"
    };

    auto fout = seqan3::sequence_file_input{std::istringstream{fastq_in}, seqan3::format_fastq{}} |
                seqan3::sequence_file_output{std::ostringstream{}, seqan3::format_fasta{}};
    fout.get_stream().flush();
    EXPECT_EQ(reinterpret_cast<std::ostringstream&>(fout.get_stream()).str(), fasta_out);
}
