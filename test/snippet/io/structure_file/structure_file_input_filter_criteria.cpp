// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges>
#include <sstream>

#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/structure_file/input.hpp>

auto input = R"(> S.cerevisiae_tRNA-PHE M10740/1-73
GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA
(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)
> example
UUGGAGUACACAACCUGUACACUCUUUC
..(((((..(((...)))..)))))... (-3.71))";

int main()
{
    seqan3::structure_file_input fin{std::istringstream{input}, seqan3::format_vienna{}};

    auto minimum_length5_filter = std::views::filter(
        [](auto const & rec)
        {
            return std::ranges::size(rec.sequence()) >= 5;
        });

    for (auto & rec : fin | minimum_length5_filter) // only record with sequence length >= 5 will "appear"
        seqan3::debug_stream << (rec.sequence() | seqan3::views::to_char) << '\n';
}
