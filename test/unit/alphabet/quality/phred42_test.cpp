// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/quality/phred42.hpp>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"
#include "phred_test_template.hpp"

using namespace seqan3;

INSTANTIATE_TYPED_TEST_CASE_P(phred42, alphabet, phred42);
INSTANTIATE_TYPED_TEST_CASE_P(phred42, alphabet_constexpr, phred42);
INSTANTIATE_TYPED_TEST_CASE_P(phred42, phred, phred42);
