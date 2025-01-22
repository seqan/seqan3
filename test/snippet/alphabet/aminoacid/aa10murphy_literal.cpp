// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/aminoacid/aa10murphy.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::aa10murphy_vector sequence1{"ACGTTA"_aa10murphy};
    seqan3::aa10murphy_vector sequence2 = "ACGTTA"_aa10murphy;
    auto sequence3 = "ACGTTA"_aa10murphy;
}
