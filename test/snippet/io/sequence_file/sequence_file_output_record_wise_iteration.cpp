// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>
#include <string>
#include <tuple>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/output.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::sequence_file_output fout{std::ostringstream{}, seqan3::format_fasta{}};

    for (int i = 0; i < 5; ++i) // ...
    {
        std::string id{"test_id"};
        seqan3::dna5_vector seq{"ACGT"_dna5};

        // ...

        fout.emplace_back(seq, id); // as individual variables
        // or:
        fout.push_back(std::tie(seq, id)); // as a tuple
    }
}
