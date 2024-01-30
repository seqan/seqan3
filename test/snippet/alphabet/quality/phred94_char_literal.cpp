// SPDX-FileCopyrightText: 2006-2023, Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2023, Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/quality/phred94.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::phred94 letter1{'!'_phred94};
    auto letter2 = '!'_phred94;
}
