// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <iostream>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

int main()
{
    using namespace seqan3::literals;

    auto char_to_rank = seqan3::to_rank('A'); // calls seqan3::custom::to_rank('A')
    static_assert(std::same_as<decltype(char_to_rank), uint8_t>);
    std::cout << static_cast<uint16_t>(char_to_rank) << '\n'; // 65

    auto dna5_to_rank = seqan3::to_rank('A'_dna5); // calls .to_char() member
    static_assert(std::same_as<decltype(dna5_to_rank), uint8_t>);
    std::cout << static_cast<uint16_t>(dna5_to_rank) << '\n'; // 0
}
