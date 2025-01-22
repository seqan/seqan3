// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <sstream>
#include <vector>

#include <seqan3/alignment/aligned_sequence/debug_stream_alignment.hpp>
#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/core/debug_stream/tuple.hpp>

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
