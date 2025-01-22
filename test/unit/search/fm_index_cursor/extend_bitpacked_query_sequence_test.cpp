// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/fm_index/fm_index.hpp>
#include <seqan3/utility/views/slice.hpp>

#include "../helper.hpp"

TEST(fm_index_cursor_test, extend_right_with_bitpacked_sequence)
{
    using namespace seqan3::literals;
    using result_t = std::vector<std::pair<uint64_t, uint64_t>>;

    seqan3::bitpacked_sequence<seqan3::dna4> seq{"ACGGTCAGGTTC"_dna4};
    seqan3::fm_index index{seq};

    auto bitpacked_query = seqan3::views::slice(seq, 1, 4);

    auto it = index.cursor();
    it.extend_right(bitpacked_query); // this line caused the compile time error
    EXPECT_EQ(seqan3::uniquify(it.locate()), (result_t{{0, 1}}));
}
