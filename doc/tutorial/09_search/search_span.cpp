// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <span>

#include <seqan3/core/debug_stream.hpp>

int main()
{
    std::string text{"Garfield the fat cat without a hat."};
    size_t start{2};
    size_t span{3};

    std::span text_view{std::data(text) + start, span}; // represent interval [2, 4]

    seqan3::debug_stream << text_view << '\n'; // Prints "rfi"
}
