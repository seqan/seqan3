// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/aliases.hpp> // includes seqan3::dna4q
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dna4_vector vec = "ACTTTGATA"_dna4;
    auto v = vec | seqan3::views::to_char;
    seqan3::debug_stream << v << '\n'; // [A,C,T,T,T,G,A,T,A]

    std::vector<seqan3::phred42> qvec{seqan3::phred42{}.assign_phred(0),
                                      seqan3::phred42{}.assign_phred(7),
                                      seqan3::phred42{}.assign_phred(5),
                                      seqan3::phred42{}.assign_phred(3),
                                      seqan3::phred42{}.assign_phred(7),
                                      seqan3::phred42{}.assign_phred(4),
                                      seqan3::phred42{}.assign_phred(30),
                                      seqan3::phred42{}.assign_phred(16),
                                      seqan3::phred42{}.assign_phred(23)};
    auto v3 = qvec | seqan3::views::to_char;
    seqan3::debug_stream << v3 << '\n'; // [!,(,&,$,(,%,?,1,8]

    std::vector<seqan3::dna4q> qcvec{{'C'_dna4, '!'_phred42},
                                     {'A'_dna4, '('_phred42},
                                     {'G'_dna4, '&'_phred42},
                                     {'T'_dna4, '$'_phred42},
                                     {'G'_dna4, '('_phred42},
                                     {'A'_dna4, '%'_phred42},
                                     {'C'_dna4, '?'_phred42},
                                     {'T'_dna4, '1'_phred42},
                                     {'A'_dna4, '8'_phred42}};
    auto v4 = qcvec | seqan3::views::to_char;
    seqan3::debug_stream << v4 << '\n'; // [C,A,G,T,G,A,C,T,A]
}
