// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>

#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/utility/views/elements.hpp>

auto input = R"(> S.cerevisiae_tRNA-PHE M10740/1-73
GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)
> example
UUGGAGUACACAACCUGUACACUCUUUC
..(((((..(((...)))..)))))... (-3.71))";

int main()
{
    using seqan3::get;

    seqan3::structure_file_input fin{std::istringstream{input},
                                     seqan3::format_vienna{},
                                     seqan3::fields<seqan3::field::id, seqan3::field::structured_seq>{}};

    // note that the order is now different, "id" comes first, because it was specified first
    for (auto & [id, struc_seq] : fin)
    {
        seqan3::debug_stream << "ID: " << id << '\n';
        // sequence and structure are part of the same vector, of type std::vector<structured_rna<rna5, wuss51>>
        // sequence and structure strings are extracted and converted to char on-the-fly
        seqan3::debug_stream << "SEQ: " << (struc_seq | seqan3::views::elements<0> | seqan3::views::to_char) << '\n';
        seqan3::debug_stream << "STRUCTURE: " << (struc_seq | seqan3::views::elements<1> | seqan3::views::to_char)
                             << '\n';
    }
}
