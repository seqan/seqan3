// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
seqan3::test::create_temporary_snippet_file my_fastq{"my.fastq",
                                                     R"//![fastq_file](
@seq1
AGCTAGCAGCGATCG
+
IIIIIHIIIIIIIII
@seq2
CGATCGATC
+
IIIIIIIII
@seq3
AGCGATCGAGGAATATAT
+
IIIIHHGIIIIHHGIIIH
)//![fastq_file]"}; // std::filesystem::current_path() / "my.fastq" will be deleted after the execution

// std::filesystem::current_path() / "output.fasta" will be deleted after the execution
seqan3::test::create_temporary_snippet_file output_fasta{"output.fasta", ""};

//![main]
#include <seqan3/io/sequence_file/all.hpp>

int main()
{
    auto current_path = std::filesystem::current_path();

    seqan3::sequence_file_output{current_path / "output.fasta"} =
        seqan3::sequence_file_input{current_path / "my.fastq"};
}
//![main]
