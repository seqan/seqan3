// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>

using seqan3::operator""_dna4;

template <typename T>
using debug_stream_test = ::testing::Test;

using alphabet_types =
    ::testing::Types<seqan3::dna4, seqan3::qualified<seqan3::dna4, seqan3::phred42>, seqan3::gapped<seqan3::dna4>>;

TYPED_TEST_SUITE(debug_stream_test, alphabet_types, );

TYPED_TEST(debug_stream_test, dna4)
{
    std::ostringstream o;
    seqan3::debug_stream_type my_stream{o};

    TypeParam val{'C'_dna4};
    my_stream << val;

    EXPECT_EQ(o.str(), "C");
}
