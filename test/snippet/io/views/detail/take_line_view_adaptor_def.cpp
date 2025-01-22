// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges>
#include <string>

#include <seqan3/core/debug_stream.hpp>

int main()
{
    //![usage]
    auto v = std::views::take_while(
        [](auto const & l)
        {
            return (l != '\r') && (l != '\n');
        });
    //![usage]
}
