// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges>
#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/range/to.hpp>
#include <seqan3/utility/views/convert.hpp>

int main()
{
    // convert from int to bool
    std::vector<int> vec{7, 5, 0, 5, 0, 0, 4, 8, -3};

    // pipe notation
    seqan3::debug_stream << (vec | seqan3::views::convert<bool>) << '\n'; // [1,1,0,1,0,0,1,1,1]

    // combinability
    seqan3::debug_stream << (vec | seqan3::views::convert<bool> | std::views::reverse) << '\n'; // [1,1,1,0,0,1,0,1,1]

    // function notation and immediate conversion to vector again
    auto bool_vec = seqan3::views::convert<bool>(vec) | seqan3::ranges::to<std::vector<bool>>();
    seqan3::debug_stream << std::boolalpha << (bool_vec == std::vector<bool>{1, 1, 0, 1, 0, 0, 1, 1, 1})
                         << '\n'; // true
}
