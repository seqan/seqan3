// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <string>
#include <vector>

#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/alphabet/views/trim_quality.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3::literals;

int main()
{
    std::vector<seqan3::phred42> vec{"II?5+"_phred42};

    // trim by phred_value
    auto v1 = vec | seqan3::views::trim_quality(20u);
    seqan3::debug_stream << v1 << '\n'; // II?5

    // trim by quality character
    auto v2 = vec | seqan3::views::trim_quality('I'_phred42);
    seqan3::debug_stream << v2 << '\n'; // II

    // function syntax
    auto v3 = seqan3::views::trim_quality(vec, '5'_phred42);
    seqan3::debug_stream << v3 << '\n'; // II?5

    // combinability
    auto v4 = seqan3::views::trim_quality(vec, 20u) | seqan3::views::to_char;
    seqan3::debug_stream << v4 << '\n'; // II?5
}
