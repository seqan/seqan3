// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/structured_rna.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/io/structure_file/output.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::structured_rna<seqan3::rna5, seqan3::wuss51> sr{'G'_rna5, '.'_wuss51};

    seqan3::structure_file_output fout{std::ostringstream{},
                                       seqan3::format_vienna{},
                                       seqan3::fields<seqan3::field::id, seqan3::field::structured_seq>{}};

    for (int i = 0; i < 10; i++) // ...
    {
        std::string id{"test_id"};
        // vector of combined data structure
        std::vector<seqan3::structured_rna<seqan3::rna5, seqan3::wuss51>> structured_sequence{sr, sr, sr, sr};

        // ...

        // note also that the order the arguments is now different, because
        // you specified that `seqan3::field::id` should be first in the fields template argument
        fout.emplace_back(id, structured_sequence);
        // or:
        fout.push_back(std::tie(id, structured_sequence));
    }
}
