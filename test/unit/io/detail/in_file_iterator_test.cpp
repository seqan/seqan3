// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/std/iterator>

#include <seqan3/io/detail/in_file_iterator.hpp>

//NOTE(h-2): This class is extensively tested via *_file_input. This is just a minimal test.

struct fake_file_t : std::vector<int>
{
    using base = std::vector<int>;

    using base::base;

    using iterator = seqan3::detail::in_file_iterator<fake_file_t>;

    size_t current_position = 0;
    bool at_end = false;

    void read_next_record()
    {
        ++current_position;
        if (current_position == size())
            at_end = true;
        else
            record_buffer = base::operator[](current_position);
    }

    iterator begin()
    {
        record_buffer = base::operator[](current_position);
        return {*this};
    }

    int record_buffer;

    fake_file_t() :
        base{}
    {}
};

TEST(in_file_iterator, concepts)
{
    using it_t = seqan3::detail::in_file_iterator<fake_file_t>;

    EXPECT_TRUE((std::input_iterator<it_t>));
}

TEST(in_file_iterator, member_types)
{
    using it_t = seqan3::detail::in_file_iterator<fake_file_t>;
    EXPECT_TRUE((std::is_same_v<typename it_t::value_type,
                                int>));
    EXPECT_TRUE((std::is_same_v<typename it_t::reference,
                                int &>));
    EXPECT_TRUE((std::is_same_v<typename it_t::const_reference,
                                int &>));
    EXPECT_TRUE((std::is_same_v<typename it_t::difference_type,
                                std::ptrdiff_t>));
    EXPECT_TRUE((std::is_same_v<typename it_t::size_type,
                                size_t>));
    EXPECT_TRUE((std::is_same_v<typename it_t::iterator_category,
                                std::input_iterator_tag>));
}

TEST(in_file_iterator, operations)
{
    using it_t = seqan3::detail::in_file_iterator<fake_file_t>;

    fake_file_t f{1, 2, 3, 4, 5, 6, 8};

    // construct
    it_t it = f.begin();

    // deref
    EXPECT_EQ(*it, 1);

    // pre-inc
    EXPECT_EQ(*(++it), 2);

    // post-inc
    it++;

    // deref
    EXPECT_EQ(*it, 3);
}

TEST(in_file_iterator, comparison)
{
    using it_t = seqan3::detail::in_file_iterator<fake_file_t>;

    fake_file_t f{1, 2, 3, 4, 5, 6, 8};
    it_t it = f.begin();

    // not at end
    EXPECT_FALSE(it == std::ranges::default_sentinel);

    // consume the entire range
    ++it; ++it; ++it; ++it; ++it; ++it; ++it;

    // at end
    EXPECT_TRUE(it == std::ranges::default_sentinel);
}
