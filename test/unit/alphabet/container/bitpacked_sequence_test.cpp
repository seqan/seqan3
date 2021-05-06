// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/composite/alphabet_variant.hpp>
#include <seqan3/alphabet/nucleotide/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/test/range/container_test_template.hpp>

INSTANTIATE_TYPED_TEST_SUITE_P(bitpacked_sequence, container_over_dna4_test, seqan3::bitpacked_sequence<seqan3::dna4>, );

using seqan3::operator""_dna4;

// https://github.com/seqan/seqan3/issues/1743
TEST(bitpacked_sequence_test, issue1743_complement_on_proxy)
{
    seqan3::bitpacked_sequence<seqan3::dna4> v{'A'_dna4};

    auto proxy = *v.begin();
    auto complement = seqan3::complement(proxy);

    EXPECT_SAME_TYPE(decltype(complement), seqan3::dna4);
    EXPECT_EQ(complement, 'T'_dna4);
}

// https://github.com/seqan/seqan3/issues/1743
TEST(bitpacked_sequence_test, issue1743_view_combinability)
{
    seqan3::bitpacked_sequence<seqan3::dna4> v{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};
    auto complement = v | seqan3::views::complement;

    EXPECT_EQ(v.size(), complement.size());
    EXPECT_RANGE_EQ(complement, (seqan3::dna4_vector{'T'_dna4, 'G'_dna4, 'C'_dna4, 'A'_dna4}));
}

// https://github.com/seqan/product_backlog/issues/371
TEST(bitpacked_sequence_test, issue371)
{
    using alphabet_t = seqan3::alphabet_variant<seqan3::dna4, seqan3::dna15>;
    seqan3::bitpacked_sequence<alphabet_t> source{};
    auto it = source.begin();
    auto end = source.end();
    it != end; // This line causes error.
}
