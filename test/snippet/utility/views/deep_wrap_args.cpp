// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/deep.hpp>

namespace my
{
inline auto const deep_take1 = seqan3::views::deep{std::views::take(1)};
}

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::dna5_vector> sequences{"AAATTT"_dna5, "CCCGGG"_dna5};

    seqan3::debug_stream << (sequences | std::views::take(1)) << '\n'; // [A,A,A,T,T,T]

    seqan3::debug_stream << (sequences | seqan3::views::deep{std::views::take(1)}) << '\n'; // [A,C]
    // constructor arguments passed via {} and arguments to underlying view hardcoded inside

    // or with an alias defined previously
    seqan3::debug_stream << (sequences | my::deep_take1) << '\n'; // [A,C]
}
