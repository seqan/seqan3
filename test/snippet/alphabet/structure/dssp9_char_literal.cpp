// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/structure/dssp9.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::dssp9 letter1{'('_dssp9};
    auto letter2 = '('_dssp9;
}
