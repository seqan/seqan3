// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <string>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/interleave.hpp>

int main()
{
    std::string u{"FOOBARBAXBAT"};
    std::string i{"in"};
    size_t s = 3;

    auto v = u | seqan3::views::interleave(s, i);

    seqan3::debug_stream << v << '\n'; // prints FOOinBARinBAXinBAT
}
