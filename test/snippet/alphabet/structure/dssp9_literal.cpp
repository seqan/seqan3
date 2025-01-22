// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/structure/dssp9.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::dssp9> sequence1{"EHHHHT"_dssp9};
    std::vector<seqan3::dssp9> sequence2 = "EHHHHT"_dssp9;
    auto sequence3 = "EHHHHT"_dssp9;
}
