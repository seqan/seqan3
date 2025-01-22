// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/validate_char_for.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    std::string_view str{"ACTTTGATAN"};
    try
    {
        seqan3::debug_stream << (str | seqan3::views::validate_char_for<seqan3::dna4>); // ACTTTGATA
    }
    catch (seqan3::invalid_char_assignment const &)
    {
        seqan3::debug_stream << "\n[ERROR] Invalid char!\n"; // Will throw on parsing 'N'
    }
}
