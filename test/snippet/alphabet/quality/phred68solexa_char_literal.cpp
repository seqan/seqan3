// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/quality/phred68solexa.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::phred68solexa letter1{'!'_phred68solexa};
    auto letter2 = '!'_phred68solexa;
}
