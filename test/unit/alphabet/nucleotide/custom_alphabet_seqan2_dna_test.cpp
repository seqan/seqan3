// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#ifdef SEQAN3_HAS_SEQAN2

// #    include "../alphabet_constexpr_test_template.hpp"
#    include "../alphabet_test_template.hpp"
// #    include "../semi_alphabet_constexpr_test_template.hpp"
#    include "../semi_alphabet_test_template.hpp"
#    include "custom_alphabet_seqan2_dna.hpp"
#    include "nucleotide_test_template.hpp"

template <>
struct nucleotide<seqan2::Dna> : public ::testing::Test
{
    static constexpr bool skip_trivial_thirdparty = true;
};

INSTANTIATE_TYPED_TEST_SUITE_P(seqan2_dna, alphabet, seqan2::Dna, );
INSTANTIATE_TYPED_TEST_SUITE_P(seqan2_dna, semi_alphabet_test, seqan2::Dna, );
// INSTANTIATE_TYPED_TEST_SUITE_P(seqan2_dna, alphabet_constexpr, seqan2::Dna, );
// INSTANTIATE_TYPED_TEST_SUITE_P(seqan2_dna, semi_alphabet_constexpr, seqan2::Dna, );
INSTANTIATE_TYPED_TEST_SUITE_P(seqan2_dna, nucleotide, seqan2::Dna, );

#endif
