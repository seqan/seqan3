// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <iostream>

#include <seqan3/core/add_enum_bitwise_operators.hpp>

enum class my_enum
{
    VAL1 = 1,
    VAL2 = 2,
    COMB = 3
};

template <>
constexpr bool seqan3::add_enum_bitwise_operators<my_enum> = true;

int main()
{
    using seqan3::operator|;

    my_enum e = my_enum::VAL1;
    my_enum e2 = e | my_enum::VAL2;

    std::cout << std::boolalpha << (e2 == my_enum::COMB) << '\n'; // true
}
