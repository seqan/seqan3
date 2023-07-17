// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <iostream>
#include <sstream>

#include <seqan3/io/stream/detail/fast_istreambuf_iterator.hpp>
#include <seqan3/io/views/detail/take_exactly_view.hpp>
#include <seqan3/io/views/detail/take_until_view.hpp>

TEST(integration_test, take_until_and_take_exactly)
{
    // This code was reduced from a failing format_sam_test
    using iterator_t = seqan3::detail::fast_istreambuf_iterator<char>;
    using sentinel_t = std::default_sentinel_t;
    using stream_view_t = std::ranges::subrange<iterator_t, sentinel_t>;

    std::istringstream stream{"HELLO WORLD"};
    std::array<char, 10> arithmetic_buffer{};

    stream_view_t stream_view{iterator_t{*stream.rdbuf()}, sentinel_t{}};
    auto stream_view_until = stream_view
                           | seqan3::detail::take_until_or_throw(
                                 [](auto)
                                 {
                                     return false;
                                 });
    auto stream_view_take2 = stream_view_until | seqan3::detail::take_exactly_or_throw(2);

    std::ranges::copy(stream_view_take2, arithmetic_buffer.data());

    EXPECT_EQ(std::string{arithmetic_buffer.data()}, "HE");
}
