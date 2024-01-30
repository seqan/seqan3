// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/aminoacid/aa10li.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::aa10li_vector sequence1{"ACGTTA"_aa10li};
    seqan3::aa10li_vector sequence2 = "ACGTTA"_aa10li;
    auto sequence3 = "ACGTTA"_aa10li;
}
