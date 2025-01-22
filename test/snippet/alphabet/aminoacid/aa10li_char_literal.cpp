// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/aminoacid/aa10li.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::aa10li letter1{'A'_aa10li};
    auto letter2 = 'A'_aa10li;
}
