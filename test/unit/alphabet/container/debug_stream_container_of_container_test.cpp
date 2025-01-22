// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/container/concatenated_sequences.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream/range.hpp>

using seqan3::operator""_dna4;

template <typename T>
class debug_stream_test : public ::testing::Test
{};

using container_of_container_types =
    ::testing::Types<std::vector<std::vector<seqan3::dna4>>,
                     seqan3::concatenated_sequences<std::vector<seqan3::dna4>>,
                     seqan3::concatenated_sequences<seqan3::bitpacked_sequence<seqan3::dna4>>>;

TYPED_TEST_SUITE(debug_stream_test, container_of_container_types, );

TYPED_TEST(debug_stream_test, container_of_container)
{
    SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY_START(-Warray-bounds, -Wstringop-overflow)
    TypeParam t1{"ACGT"_dna4, "ACGT"_dna4, "GAGGA"_dna4};
    SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY_STOP

    std::ostringstream o;
    seqan3::debug_stream_type my_stream{o};

    my_stream << TypeParam{};

    o.flush();
    EXPECT_EQ(o.str(), "[]");

    my_stream << ", " << t1;

    o.flush();
    EXPECT_EQ(o.str(), "[], [ACGT,ACGT,GAGGA]");
}
