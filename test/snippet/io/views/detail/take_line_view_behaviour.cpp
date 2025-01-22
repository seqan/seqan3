// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <string>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/views/detail/take_line_view.hpp>

int main()
{
    std::string vec{"foo\nbar"};
    auto v = vec | seqan3::detail::take_line;
    seqan3::debug_stream << v << '\n'; // [f,o,o]

    auto v2 = vec | std::views::reverse | seqan3::detail::take_line;
    seqan3::debug_stream << v2 << '\n'; // [r,a,b]
    seqan3::debug_stream << v2 << '\n'; // [r,a,b] (parsing it again gives us the same result)
}
