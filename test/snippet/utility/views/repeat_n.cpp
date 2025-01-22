// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges>
#include <string>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/repeat_n.hpp>

int main()
{
    auto v = seqan3::views::repeat_n(std::string{"foo"}, 5);

    seqan3::debug_stream << v.size() << '\n'; // prints 5
    seqan3::debug_stream << v << '\n';        // prints ["foo", "foo", "foo", "foo", "foo"]

    v[0] = std::string{"foobar"};

    // Now it always prints "foobar"
    seqan3::debug_stream << *std::ranges::begin(v) << '\n'; // prints "foobar"
    seqan3::debug_stream << v.size() << '\n';               // prints 5; Note: the size cannot be changed anymore
}
