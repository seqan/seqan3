// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

// Checks that defining simd_dna4 works without putting it into the seqan3 namespace.

#define SEQAN3_USE_NAMESPACE 0
#include <seqan3/test/performance/simd_dna4.hpp>

int main()
{
    simd_dna4 letter{};
}
