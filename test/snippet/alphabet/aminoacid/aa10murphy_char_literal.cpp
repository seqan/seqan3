// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/aminoacid/aa10murphy.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::aa10murphy letter1{'A'_aa10murphy};
    auto letter2 = 'A'_aa10murphy;
}
