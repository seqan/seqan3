// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
// std::filesystem::current_path() / "my.fasta" will be deleted after the execution
seqan3::test::create_temporary_snippet_file my_fasta{"my.fasta", ""};

//![main]
#include <filesystem>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>

int main()
{
    using namespace seqan3::literals;

    auto fasta_file = std::filesystem::current_path() / "my.fasta";

    {
        // Create a ./my.fasta file.
        seqan3::sequence_file_output fout{fasta_file};
        fout.emplace_back("ACGT"_dna4, "TEST1");
        fout.emplace_back("AGGCTGA"_dna4, "Test2");
        fout.emplace_back("GGAGTATAATATATATATATATAT"_dna4, "Test3");
    }

    // FASTA with DNA sequences assumed, regular std::ifstream taken as stream
    seqan3::sequence_file_input fin{fasta_file};
}
//![main]
