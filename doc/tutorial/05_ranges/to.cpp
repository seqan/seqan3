// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <deque>
#include <forward_list>
#include <list>
#include <vector>

#include <seqan3/utility/range/to.hpp>

int main()
{
    auto lst = std::views::iota(1, 10); // some range over the numbers 1-10

    // convert range to vector using pipe syntax
    auto vec0 = lst | seqan3::ranges::to<std::vector<int>>();
    static_assert(std::same_as<decltype(vec0), std::vector<int>>);

    // convert range to vector but auto deducing the element type
    auto vec1 = lst | seqan3::ranges::to<std::vector>();
    static_assert(std::same_as<decltype(vec1), std::vector<int>>);

    // convert range to vector using function call syntax
    auto vec2 = seqan3::ranges::to<std::vector<int>>(lst);
    static_assert(std::same_as<decltype(vec2), std::vector<int>>);

    // using function call syntax and auto deducing element type
    auto vec3 = seqan3::ranges::to<std::vector>(lst);
    static_assert(std::same_as<decltype(vec3), std::vector<int>>);

    // convert nested ranges into nested containers
    auto nested_lst = std::list<std::forward_list<int>>{{1, 2, 3}, {4, 5, 6, 7}};
    auto vec4 = nested_lst | seqan3::ranges::to<std::vector<std::vector<int>>>();
    static_assert(std::same_as<decltype(vec4), std::vector<std::vector<int>>>);

    // different supported container types
    auto vec5 = lst | seqan3::ranges::to<std::list>();
    static_assert(std::same_as<decltype(vec5), std::list<int>>);

    auto vec6 = lst | seqan3::ranges::to<std::deque>();
    static_assert(std::same_as<decltype(vec6), std::deque<int>>);
}
