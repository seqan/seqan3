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
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

template <typename T>
using alphabet_debug_stream = ::testing::Test;

using test_types = ::testing::Types<dna4, qualified<dna4, phred42>, gapped<dna4>>;

TYPED_TEST_CASE(alphabet_debug_stream, test_types);

TYPED_TEST(alphabet_debug_stream, debug_streaming)
{
    std::ostringstream o;
    debug_stream_type my_stream{o};

    TypeParam val{'C'_dna4};
    my_stream << val;

    o.flush();
    EXPECT_EQ(o.str().size(), 1u);
    EXPECT_EQ(to_char(o.str()[0]), 'C');
}
