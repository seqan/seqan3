// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

//![main]
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

int main()
{
    seqan3::sequence_file_input{"input.fastq"} | seqan3::sequence_file_output{"output.fasta"};
    // no variables are created here, all input is immediately streamed to the output
    // the format is converted on-the-fly
}
//![main]

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
// std::filesystem::current_path() / "input.fastq" will be deleted after the execution
seqan3::test::create_temporary_snippet_file input_fastq{"input.fastq", "\n"};
// std::filesystem::current_path() / "output.fasta" will be deleted after the execution
seqan3::test::create_temporary_snippet_file output_fasta{"output.fasta", ""};
