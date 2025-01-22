// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>
#include <string>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::sequence_file_output fout{std::ostringstream{}, seqan3::format_fasta{}};
    for (int i = 0; i < 5; ++i) // some criteria
    {
        seqan3::record<seqan3::type_list<seqan3::dna5_vector, std::string>,
                       seqan3::fields<seqan3::field::seq, seqan3::field::id>>
            r{"ACGT"_dna5, "ID1"};

        // ...

        fout.push_back(r);
    }
}
