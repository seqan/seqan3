// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/structure/all.hpp>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"

using namespace seqan3;

using structured_rna_types = ::testing::Types<structured_rna<rna5, dot_bracket3>, structured_rna<rna4, wuss51>>;

INSTANTIATE_TYPED_TEST_CASE_P(structured_rna, alphabet, structured_rna_types);
INSTANTIATE_TYPED_TEST_CASE_P(structured_rna, alphabet_constexpr, structured_rna_types);
