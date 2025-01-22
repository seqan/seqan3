// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

//![main]
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>

int main()
{
    seqan3::sequence_file_input fin{"my.fasta"};

    // iterate over every record in the file and split it into fields:
    for (auto & [seq, id, qual] : fin)
    {
        // print the fields:
        seqan3::debug_stream << "ID:  " << id << '\n';
        seqan3::debug_stream << "SEQ: " << seq << '\n';
        seqan3::debug_stream << "QUAL:" << qual << '\n'; // qual is empty for FASTA files
    }
}
//![main]

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
// std::filesystem::current_path() / "my.fasta" will be deleted after the execution
seqan3::test::create_temporary_snippet_file my_fastq{"my.fasta", "\n"};
