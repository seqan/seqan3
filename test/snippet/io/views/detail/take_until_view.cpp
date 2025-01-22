// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges> // for std::views::reverse
#include <string>

#include <seqan3/core/debug_stream.hpp>                 // for debug_stream
#include <seqan3/io/views/detail/take_until_view.hpp>   // for detail::take_until*
#include <seqan3/utility/char_operations/predicate.hpp> // for is_char
#include <seqan3/utility/views/single_pass_input.hpp>   // for views::single_pass_input

int main()
{
    // regular usage
    std::string vec{"foo\nbar"};
    auto v = vec | seqan3::detail::take_until(seqan3::is_char<'\n'>); // or use a lambda
    seqan3::debug_stream << v << '\n';                                // "foo"

    auto v2 = vec | std::views::reverse | seqan3::detail::take_until(seqan3::is_char<'\n'>);
    seqan3::debug_stream << v2 << '\n'; // "rab"

    // consuming behaviour
    std::string vec2{"foo      bar"}; // ← multiple spaces
    auto vin = vec2 | seqan3::views::single_pass_input;
    auto v3 = vin | seqan3::detail::take_until_and_consume(seqan3::is_blank);
    seqan3::debug_stream << v3 << '\n';                       // "foo"
    seqan3::debug_stream << *std::ranges::begin(vin) << '\n'; // "b", the spaces where skipped
}
