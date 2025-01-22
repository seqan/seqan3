// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/mask/mask.hpp>

TEST(debug_stream_test, mask)
{
    {
        std::ostringstream o;
        seqan3::debug_stream_type my_stream{o};
        my_stream << seqan3::mask::masked;
        EXPECT_EQ(o.str(), "MASKED");
    }

    {
        std::ostringstream o;
        seqan3::debug_stream_type my_stream{o};
        my_stream << seqan3::mask::unmasked;
        EXPECT_EQ(o.str(), "UNMASKED");
    }
}
