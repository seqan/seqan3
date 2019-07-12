// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/test/cereal.hpp>

using namespace seqan3;

template <typename T>
using alphabet_cereal = ::testing::Test;

using test_types = ::testing::Types<dna4, qualified<dna4, phred42>, gapped<dna4>>;

TYPED_TEST_CASE(alphabet_cereal, test_types);

TYPED_TEST(alphabet_cereal, serialisation)
{
    TypeParam letter;

    assign_rank_to(1 % alphabet_size<TypeParam>, letter);
    test::do_serialisation(letter);

    std::vector<TypeParam> vec;
    vec.resize(10);
    for (unsigned i = 0; i < 10; ++i)
        assign_rank_to(i % alphabet_size<TypeParam>, vec[i]);
    test::do_serialisation(vec);
}
