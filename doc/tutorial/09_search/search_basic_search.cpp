// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/core/debug_stream.hpp> // pretty printing
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/search/search.hpp>

using namespace std::string_literals; // for using the ""s string literal

int main()
{
    std::string text{"Garfield the fat cat without a hat."};
    seqan3::fm_index index{text};
    seqan3::debug_stream << search("cat"s, index) << '\n'; // [<query_id:0, reference_id:0, reference_pos:17>]
}
