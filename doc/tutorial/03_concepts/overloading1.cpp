// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <concepts>
#include <iostream> // for std::cout

template <std::integral t>
void print(t const v)
{
    std::cout << "integral value: " << v << '\n';
}

int main()
{
    int i{4};
    unsigned u{3};

    print(i); // prints "integral value: 4"
    print(u); // prints "integral value: 3"
}
