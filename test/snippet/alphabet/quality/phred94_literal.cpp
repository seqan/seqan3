// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/alphabet/quality/phred94.hpp>

int main()
{
    using namespace seqan3::literals;

    std::vector<seqan3::phred94> sequence1{"##!!##"_phred94};
    std::vector<seqan3::phred94> sequence2 = "##!!##"_phred94;
    auto sequence3 = "##!!##"_phred94;
}
