// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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
