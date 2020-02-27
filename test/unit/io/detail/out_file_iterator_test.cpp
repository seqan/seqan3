// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/std/iterator>

#include <seqan3/io/detail/out_file_iterator.hpp>

//NOTE(h-2): This class is extensively tested via *_file_output. This is just a minimal test.

TEST(out_file_iterator, concepts)
{
    using it_t = seqan3::detail::out_file_iterator<std::vector<int>>;

    EXPECT_TRUE((std::output_iterator<it_t, int>));
}

TEST(out_file_iterator, member_types)
{
    using it_t = seqan3::detail::out_file_iterator<std::vector<int>>;
    EXPECT_TRUE((std::is_same_v<typename it_t::value_type,
                                void>));
    EXPECT_TRUE((std::is_same_v<typename it_t::reference,
                                void>));
    EXPECT_TRUE((std::is_same_v<typename it_t::const_reference,
                                void>));
    EXPECT_TRUE((std::is_same_v<typename it_t::difference_type,
                                std::ptrdiff_t>));
    EXPECT_TRUE((std::is_same_v<typename it_t::size_type,
                                void>));
    EXPECT_TRUE((std::is_same_v<typename it_t::iterator_category,
                                std::output_iterator_tag>));
}

TEST(out_file_iterator, operations)
{
    using it_t = seqan3::detail::out_file_iterator<std::vector<int>>;

    std::vector<int> fake_file;

    // construct
    it_t it{fake_file};

    // pre-inc, no-op
    ++it;

    // post-inc, no-op
    it++;

    // assign to iterator
    it = 3;
    EXPECT_EQ(fake_file.size(), 1u);
    EXPECT_EQ(fake_file[0], 3);

    // assign to deref'ed iterator
    *it = 7;
    EXPECT_EQ(fake_file.size(), 2u);
    EXPECT_EQ(fake_file[0], 3);
    EXPECT_EQ(fake_file[1], 7);

    // assign to deref'ed iterator-increment
    *it++ = 9;
    EXPECT_EQ(fake_file.size(), 3u);
    EXPECT_EQ(fake_file[0], 3);
    EXPECT_EQ(fake_file[1], 7);
    EXPECT_EQ(fake_file[2], 9);
}

TEST(out_file_iterator, comparison)
{
    using it_t = seqan3::detail::out_file_iterator<std::vector<int>>;

    std::vector<int> fake_file;
    it_t it{fake_file};

    // never at end
    EXPECT_FALSE(it == std::ranges::default_sentinel);
}
