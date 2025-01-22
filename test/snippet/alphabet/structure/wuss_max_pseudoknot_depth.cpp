// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <iostream>

#include <seqan3/alphabet/structure/wuss.hpp>

int main()
{
    uint8_t max_depth_member = seqan3::wuss51::max_pseudoknot_depth;
    uint8_t max_depth_meta = seqan3::max_pseudoknot_depth<seqan3::wuss51>;
    std::cout << static_cast<uint16_t>(max_depth_member) << '\n'; // 22
    std::cout << static_cast<uint16_t>(max_depth_meta) << '\n';   // 22
}
