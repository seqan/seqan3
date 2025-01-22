// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <ranges>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/deep.hpp>

namespace my
{
// You can create a permanent alias:
inline auto const deep_reverse = seqan3::views::deep{std::views::reverse};
} // namespace my

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::dna5_vector> sequences{"AAATTT"_dna5, "CCCGGG"_dna5};

    seqan3::debug_stream << (sequences | std::views::reverse) << '\n'; // [CCCGGG,AAATTT]

    // These two are equivalent:
    seqan3::debug_stream << (sequences | my::deep_reverse) << '\n';                         // [TTTAAA,GGGCCC]
    seqan3::debug_stream << (sequences | seqan3::views::deep{std::views::reverse}) << '\n'; // [TTTAAA,GGGCCC]
}
