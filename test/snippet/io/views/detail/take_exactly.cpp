// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <string>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/views/detail/take_exactly_view.hpp>

int main()
{
    std::string vec{"foobar"};
    auto v = vec | seqan3::detail::take_exactly(3);       // or seqan3::detail::take_exactly_or_throw
    seqan3::debug_stream << v << '\n';                    // "foo"
    seqan3::debug_stream << std::ranges::size(v) << '\n'; // 3

    auto v2 = vec | seqan3::detail::take_exactly(9);
    seqan3::debug_stream << std::ranges::size(v2) << '\n'; // 9 <- here be dragons! (undefined behaviour)
}
