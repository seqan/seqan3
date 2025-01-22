// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <iostream>
#include <ranges> // include all of the standard library's views
#include <vector>

int main()
{
    std::vector vec{1, 2, 3, 4, 5, 6};
    auto v = vec
           | std::views::filter(
                 [](auto const i)
                 {
                     return i % 2 == 0;
                 })
           | std::views::transform(
                 [](auto const i)
                 {
                     return i * i;
                 });

    std::cout << *v.begin() << '\n'; // prints 4
}
