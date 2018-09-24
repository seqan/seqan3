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

#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/view/take.hpp>
#include <seqan3/test/tmp_filename.hpp>
#include <seqan3/std/iterator>

using namespace seqan3;
using namespace seqan3::literal;

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

    sequence_file_input fin{std::istringstream{input}, sequence_file_format_fasta{}};
    sequence_file_output fout{std::ostringstream{}, sequence_file_format_fasta{}};
    fout.options.fasta_letters_per_line = 0;

    fout = fin;

    fout.get_stream().flush();
    EXPECT_EQ(fout.get_stream().str(), output_comp);
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
    sequence_file_input{std::istringstream{input}, sequence_file_format_fasta{}} |
        sequence_file_output{std::ostringstream{}, sequence_file_format_fasta{}};

    // valid with assignment and check contents
    auto fout = sequence_file_input{std::istringstream{input}, sequence_file_format_fasta{}} |
                sequence_file_output{std::ostringstream{}, sequence_file_format_fasta{}};

    fout.get_stream().flush();
    EXPECT_EQ(fout.get_stream().str(), input);
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
    sequence_file_input{std::istringstream{input}, sequence_file_format_fasta{}} | view::persist | view::take(2) |
        sequence_file_output{std::ostringstream{}, sequence_file_format_fasta{}};

    // valid with assignment and check contents
    auto fout = sequence_file_input{std::istringstream{input}, sequence_file_format_fasta{}}
              | view::persist
              | view::take(2)
              | sequence_file_output{std::ostringstream{}, sequence_file_format_fasta{}};

    fout.get_stream().flush();
    EXPECT_EQ(fout.get_stream().str(), output);
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

    auto fout = sequence_file_input{std::istringstream{fastq_in}, sequence_file_format_fastq{}} |
                sequence_file_output{std::ostringstream{}, sequence_file_format_fasta{}};
    fout.get_stream().flush();
    EXPECT_EQ(fout.get_stream().str(), fasta_out);
}
