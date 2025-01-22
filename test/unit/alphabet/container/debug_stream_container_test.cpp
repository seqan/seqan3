// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream/range.hpp>
#include <seqan3/utility/container/small_vector.hpp>

using seqan3::operator""_dna4;

template <typename T>
class debug_stream_test : public ::testing::Test
{};

using container_types = ::testing::Types<std::vector<seqan3::dna4>,
                                         seqan3::bitpacked_sequence<seqan3::dna4>,
                                         seqan3::small_vector<seqan3::dna4, 1000>>;

TYPED_TEST_SUITE(debug_stream_test, container_types, );

TYPED_TEST(debug_stream_test, container)
{
    TypeParam t1{'A'_dna4, 'C'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4};

    std::ostringstream o;
    seqan3::debug_stream_type my_stream{o};

    my_stream << TypeParam{};

    o.flush();
    EXPECT_EQ(o.str(), "");

    my_stream << t1;

    o.flush();
    EXPECT_EQ(o.str(), "ACCGT");
}
