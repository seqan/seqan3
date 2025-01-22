// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/pairwise_combine.hpp>

int main()
{
    std::vector vec{'a', 'b', 'c', 'a'};
    for (auto res : vec | seqan3::views::pairwise_combine)
    {
        seqan3::debug_stream << res << '\n';
    }

    // Possible Output:
    // (a,b)
    // (a,c)
    // (a,a)
    // (b,c)
    // (b,a)
    // (c,a)
}
