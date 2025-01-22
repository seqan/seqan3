// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <iostream>

#include <seqan3/alphabet/structure/wuss.hpp>

int main()
{
    using namespace seqan3::literals;

    bool is_unpaired_char_member = '.'_wuss51.is_unpaired();
    bool is_unpaired_char_free = seqan3::is_unpaired('{'_wuss51);

    std::cout << std::boolalpha << is_unpaired_char_member << '\n'; // true
    std::cout << std::boolalpha << is_unpaired_char_free << '\n';   // false
}
