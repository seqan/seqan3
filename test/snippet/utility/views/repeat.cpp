// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/repeat.hpp>

int main()
{
    auto v = seqan3::views::repeat('A');

    seqan3::debug_stream << *std::ranges::begin(v) << '\n'; // prints 'A'
    seqan3::debug_stream << v[12355] << '\n';               // also prints 'A'. It always prints 'A'

    v[1345] = 'C';

    // Now it always prints 'C'
    seqan3::debug_stream << *std::ranges::begin(v) << '\n'; // prints 'C'
    seqan3::debug_stream << v[12355] << '\n';               // prints 'C'
}
