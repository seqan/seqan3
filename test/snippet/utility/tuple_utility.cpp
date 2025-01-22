// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/tuple/split.hpp>

int main()
{
    // Split at position 2.
    std::tuple<int, char, float, std::string> t{1, 'c', 0.3, "hello"};
    auto [left, right] = seqan3::tuple_split<2>(t);
    static_assert(std::same_as<decltype(left), std::tuple<int, char>>);
    static_assert(std::same_as<decltype(right), std::tuple<float, std::string>>);

    // Split at position 0.
    auto [left1, right1] = seqan3::tuple_split<0>(t);
    static_assert(std::same_as<decltype(left1), std::tuple<>>);
    static_assert(std::same_as<decltype(right1), std::tuple<int, char, float, std::string>>);

    // Split at position 4.
    auto [left2, right2] = seqan3::tuple_split<4>(t);
    static_assert(std::same_as<decltype(left2), std::tuple<int, char, float, std::string>>);
    static_assert(std::same_as<decltype(right2), std::tuple<>>);
}
