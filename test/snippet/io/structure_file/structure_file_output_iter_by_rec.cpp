// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/io/structure_file/output.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::structure_file_output fout{std::ostringstream{}, seqan3::format_vienna{}};

    for (int i = 0; i < 10; i++) // ...
    {
        std::string id{"test_id"};
        seqan3::rna5_vector seq{"ACGU"_rna5};
        std::vector<seqan3::wuss51> structure{".()."_wuss51};

        // ...

        fout.emplace_back(seq, id, structure); // as individual variables
        // or:
        fout.push_back(std::tie(seq, id, structure)); // as a tuple
    }
}
