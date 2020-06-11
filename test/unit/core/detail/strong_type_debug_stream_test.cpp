// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
