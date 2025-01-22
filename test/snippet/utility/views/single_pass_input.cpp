// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <string>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/single_pass_input.hpp>

int main()
{
    std::string str{"hello"};
    auto v = str | seqan3::views::single_pass_input;
    auto b = v.begin();

    seqan3::debug_stream << *b << '\n';     // prints 'h'
    seqan3::debug_stream << *(++b) << '\n'; // prints 'e'
    seqan3::debug_stream << *(++b) << '\n'; // prints 'l'
    seqan3::debug_stream << *(++b) << '\n'; // prints 'l'
    seqan3::debug_stream << *(++b) << '\n'; // prints 'o'
}
