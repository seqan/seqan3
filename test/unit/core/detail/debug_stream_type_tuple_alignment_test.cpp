// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <sstream>
#include <vector>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alignment/aligned_sequence/debug_stream_alignment.hpp>
#include <seqan3/core/detail/debug_stream_tuple.hpp>

// https://github.com/seqan/product_backlog/issues/125
TEST(debug_stream_type, issue_125)
{
    using aligned_sequence_t = std::vector<seqan3::gapped<char>>;
    using alignment_t = std::pair<aligned_sequence_t, aligned_sequence_t>;

    std::ostringstream oss;
    seqan3::debug_stream_type stream{oss};
    stream << alignment_t{};
    EXPECT_EQ(oss.str(), "");
}
