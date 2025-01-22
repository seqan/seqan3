// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/core/detail/strong_type.hpp>

struct my_type : seqan3::detail::strong_type<int, my_type>
{
    using seqan3::detail::strong_type<int, my_type>::strong_type;
};

TEST(strong_type, debug_stremable)
{
    my_type obj{10};

    std::ostringstream buffer{};
    seqan3::debug_stream_type stream{buffer};
    stream << obj;
    EXPECT_EQ(buffer.str(), "10");
}
