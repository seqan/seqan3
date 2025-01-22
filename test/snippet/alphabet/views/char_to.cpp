// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    std::string str{"ACTTTGATAN"};
    seqan3::debug_stream << (str | seqan3::views::char_to<seqan3::dna4>) << '\n'; // ACTTTGATAA
    seqan3::debug_stream << (str | seqan3::views::char_to<seqan3::dna5>) << '\n'; // ACTTTGATAN
}
