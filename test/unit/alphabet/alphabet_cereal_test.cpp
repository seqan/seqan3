// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/test/cereal.hpp>

template <typename t>
using alphabet_cereal = ::testing::Test;

using test_types =
    ::testing::Types<seqan3::dna4, seqan3::qualified<seqan3::dna4, seqan3::phred42>, seqan3::gapped<seqan3::dna4>>;

TYPED_TEST_SUITE(alphabet_cereal, test_types, );

TYPED_TEST(alphabet_cereal, serialisation)
{
    TypeParam letter;

    seqan3::assign_rank_to(1 % seqan3::alphabet_size<TypeParam>, letter);
    seqan3::test::do_serialisation(letter);

    std::vector<TypeParam> vec;
    vec.resize(10);
    for (unsigned i = 0; i < 10; ++i)
        seqan3::assign_rank_to(i % seqan3::alphabet_size<TypeParam>, vec[i]);
    seqan3::test::do_serialisation(vec);
}
