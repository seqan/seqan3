// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <filesystem>
#include <sstream>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

auto input = R"(@TEST1
ACGT
+
##!#
@Test2
AGGCTGA
+
##!#!!!
@Test3
GGAGTATAATATATATATATATAT
+
##!###!###!###!###!###!#)";

int main()
{
    // file format conversion in one line:
    seqan3::sequence_file_output{std::ostringstream{}, seqan3::format_fasta{}} =
        seqan3::sequence_file_input{std::istringstream{input}, seqan3::format_fastq{}};

    // with seqan3::sequence_file_output as a variable:
    seqan3::sequence_file_output fout{std::ostringstream{}, seqan3::format_fasta{}};
    seqan3::sequence_file_input fin{std::istringstream{input}, seqan3::format_fastq{}};
    fout = fin;

    // or in pipe notation:
    seqan3::sequence_file_input{std::istringstream{input}, seqan3::format_fastq{}}
        | seqan3::sequence_file_output{std::ostringstream{}, seqan3::format_fasta{}};
}
