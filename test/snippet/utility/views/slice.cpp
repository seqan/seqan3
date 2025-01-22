// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges> // provides std::views::reverse
#include <string>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/slice.hpp> // provides views::slice

int main()
{
    std::string s{"foobar"};
    auto v = s | seqan3::views::slice(1, 4);
    seqan3::debug_stream << v << '\n'; // "oob"

    auto v2 = s | std::views::reverse | seqan3::views::slice(1, 4);
    seqan3::debug_stream << v2 << '\n'; // "abo"
}
