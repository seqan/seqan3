// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/cigar/cigar.hpp>

TEST(debug_stream_test, cigar)
{
    std::ostringstream o;
    seqan3::debug_stream_type my_stream{o};

    seqan3::cigar c1{};
    c1.assign_string("223M");

    my_stream << c1;

    o.flush();
    EXPECT_EQ(o.str().size(), 4u);
    EXPECT_EQ(o.str(), "223M");
}
