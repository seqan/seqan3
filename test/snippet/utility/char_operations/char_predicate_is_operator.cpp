// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <iostream>

#include <seqan3/utility/char_operations/predicate.hpp>

int main()
{
    char chr{'1'};
    constexpr auto my_cond = seqan3::is_char<'%'> || seqan3::is_digit;
    bool is_percent = my_cond(chr);
    std::cout << std::boolalpha << is_percent << '\n'; // true
}
