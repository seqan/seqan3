// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>

#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/io/structure_file/output.hpp>

auto input = R"(> S.cerevisiae_tRNA-PHE M10740/1-73
GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)
> example
UUGGAGUACACAACCUGUACACUCUUUC
..(((((..(((...)))..)))))... (-3.71))";

int main()
{
    bool criteria = true;
    seqan3::structure_file_input fin{std::istringstream{input},
                                     seqan3::format_vienna{},
                                     seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::structure>{}};
    // the output doesn't have to match the configuration of the input
    seqan3::structure_file_output fout{std::ostringstream{}, seqan3::format_vienna{}};

    for (auto & r : fin)
    {
        if (criteria) // r fulfills some filter criterium
            fout.push_back(r);
    }
}
