// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
// std::filesystem::current_path() / "my.fasta" will be deleted after the execution
seqan3::test::create_temporary_snippet_file my_fasta{"my.fasta", ""};

//![main]
#include <filesystem>

#include <seqan3/io/sequence_file/output.hpp>

int main()
{
    auto fasta_file = std::filesystem::current_path() / "my.fasta";

    // FASTA format detected, std::ofstream opened for file
    seqan3::sequence_file_output fin{fasta_file};
}
//![main]
