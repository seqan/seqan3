// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

//![start]
#include <iostream>
#include <ranges>

#include <seqan3/alphabet/nucleotide/all.hpp>

using seqan3::operator""_dna5;

//![start]
auto const my_convert_to_char_view = std::views::transform(
    [](auto const alph)
    {
        return seqan3::to_char(alph);
    });

//![end]
int main()
{
    std::vector<seqan3::dna5> vec{"ATTAGATTA"_dna5};
    // std::cout << vec[0] << '\n';                 // won't work

    auto v = vec | my_convert_to_char_view;

    std::cout << v[0] << '\n'; // prints "A"
}
//![end]
