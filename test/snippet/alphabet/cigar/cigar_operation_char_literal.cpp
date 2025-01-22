// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/cigar/cigar.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::cigar::operation letter1{'M'_cigar_operation};
    auto letter2 = 'M'_cigar_operation;
}
