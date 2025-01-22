// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/nucleotide/dna4.hpp>

int main()
{
    // does not work:
    // seqan3::dna4 my_letter{0};      // we want to set the default, an A
    // seqan3::dna4 my_letter{'A'};    // we also want to set an A, but we are setting value 65

    // std::cout << my_letter;         // you expect 'A', but how would you access the number?
}
