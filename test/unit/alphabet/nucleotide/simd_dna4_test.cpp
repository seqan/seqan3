// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/core/debug_stream/range.hpp>
#include <seqan3/test/performance/simd_dna4.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "nucleotide_test_template.hpp"

INSTANTIATE_TYPED_TEST_SUITE_P(simdifyable_dna4, alphabet, seqan3::simd_dna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(simdifyable_dna4, semi_alphabet_test, seqan3::simd_dna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(simdifyable_dna4, alphabet_constexpr, seqan3::simd_dna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(simdifyable_dna4, semi_alphabet_constexpr, seqan3::simd_dna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(simdifyable_dna4, nucleotide, seqan3::simd_dna4, );

TEST(simdifyable_dna4, to_char_assign_char)
{
    EXPECT_EQ(seqan3::to_char(seqan3::simd_dna4{}.assign_char('A')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::simd_dna4{}.assign_char('C')), 'C');
    EXPECT_EQ(seqan3::to_char(seqan3::simd_dna4{}.assign_char('G')), 'G');
    EXPECT_EQ(seqan3::to_char(seqan3::simd_dna4{}.assign_char('T')), 'T');

    EXPECT_EQ(seqan3::to_char(seqan3::simd_dna4{}.assign_char('a')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::simd_dna4{}.assign_char('c')), 'C');
    EXPECT_EQ(seqan3::to_char(seqan3::simd_dna4{}.assign_char('g')), 'G');
    EXPECT_EQ(seqan3::to_char(seqan3::simd_dna4{}.assign_char('t')), 'T');
}
