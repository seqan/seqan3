// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/alphabet/mask/masked.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

using masked_types = ::testing::Types<seqan3::masked<seqan3::dna4>, seqan3::masked<seqan3::dna5>>;

INSTANTIATE_TYPED_TEST_SUITE_P(masked, alphabet, masked_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(masked, semi_alphabet_test, masked_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(masked, alphabet_constexpr, masked_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(masked, semi_alphabet_constexpr, masked_types, );
