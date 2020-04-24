// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/mask/masked.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/aminoacid/aa20.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

using masked_types = ::testing::Types<seqan3::masked<seqan3::dna4>, seqan3::masked<seqan3::dna5>>;

INSTANTIATE_TYPED_TEST_SUITE_P(masked, alphabet_, masked_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(masked, semi_alphabet_test, masked_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(masked, alphabet_constexpr, masked_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(masked, semi_alphabet_constexpr, masked_types, );
