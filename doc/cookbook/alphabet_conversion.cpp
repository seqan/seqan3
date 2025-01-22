// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/core/debug_stream.hpp>

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;
using seqan3::operator""_phred42;

int main()
{
    // A vector of combined sequence and quality information.
    std::vector<seqan3::dna4q> sequence1{{'A'_dna4, '!'_phred42},
                                         {'C'_dna4, 'A'_phred42},
                                         {'G'_dna4, '6'_phred42},
                                         {'T'_dna4, '&'_phred42}};
    // A vector of dna5.
    std::vector<seqan3::dna5> sequence2{"AGNCGTNNCAN"_dna5};

    // Convert dna4q to dna4.
    // Since `sequence1` is an lvalue, we capture `in` via const &. When unsure, use the general case below.
    auto view1 = sequence1
               | std::views::transform(
                     [](auto const & in)
                     {
                         return static_cast<seqan3::dna4>(in);
                     });
    seqan3::debug_stream << view1 << '\n'; // ACGT

    // Convert dna5 to dna4.
    // General case: Perfect forward.
    auto view2 = sequence2 | std::views::take(8)
               | std::views::transform(
                     [](auto && in)
                     {
                         return static_cast<seqan3::dna4>(std::forward<decltype(in)>(in));
                     });
    seqan3::debug_stream << view2 << '\n'; // AGACGTAA

    return 0;
}
