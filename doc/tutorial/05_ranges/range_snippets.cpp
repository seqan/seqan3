// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <iostream>
#include <ranges> // include all of the standard library's views
#include <vector>

int main()
{
    {
        //![all]
        //![def]
        std::vector vec{1, 2, 3, 4, 5, 6};
        //![rev_def]
        auto v = std::views::reverse(vec);
        //![rev_def]
        //![def]

        std::cout << *v.begin() << '\n';
        //![all]
    }

    {
        //![assign_through]
        //![piped]
        std::vector vec{1, 2, 3, 4, 5, 6};
        auto v = vec | std::views::reverse | std::views::drop(2);

        std::cout << *v.begin() << '\n';
        //![piped]
        *v.begin() = 42; // now vec == {1, 2, 3, 42, 5, 6 } !!
        //![assign_through]
    }
}
