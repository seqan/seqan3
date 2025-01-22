// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/structure/wuss.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::wuss51> sequence1{".<..>."_wuss51};
    std::vector<seqan3::wuss51> sequence2 = ".<..>."_wuss51;
    auto sequence3 = ".<..>."_wuss51;
}
