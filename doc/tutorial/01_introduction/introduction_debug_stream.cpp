// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

//! [debug]
#include <iostream> // for std::cerr
#include <vector>   // for std::vector

#include <seqan3/core/debug_stream.hpp> // for debug_stream

int main()
{
    std::vector<int> vec{-1, 0, 1};
    seqan3::debug_stream << vec << '\n'; // => [-1,0,1]
    // std::cerr << vec << '\n';            // compiler error: no operator<< for std::vector<int>
}
//! [debug]
