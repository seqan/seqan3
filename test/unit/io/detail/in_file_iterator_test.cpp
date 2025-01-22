// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iostream>
#include <iterator>
#include <memory>

#include <seqan3/io/detail/in_file_iterator.hpp>
#include <seqan3/io/detail/misc_input.hpp>
#include <seqan3/test/expect_same_type.hpp>

struct fake_file_t
{
    using value_type = char;
    using reference = char &;
    using size_type = size_t;
    using difference_type = std::ptrdiff_t;

    using iterator = seqan3::detail::in_file_iterator<fake_file_t>;

    using stream_ptr_t = std::unique_ptr<std::basic_istream<char>>;

    bool at_end{false};
    value_type record_buffer{};
    stream_ptr_t secondary_stream{};
    std::streampos position_buffer{};

    void read_next_record()
    {
        position_buffer = secondary_stream->tellg();

        // at end if we could not read further
        if (secondary_stream->eof())
        {
            at_end = true;
            return;
        }

        secondary_stream->get(record_buffer);
    }

    iterator begin()
    {
        this->read_next_record();
        return {*this};
    }

    fake_file_t(std::istringstream in) : secondary_stream{new std::istringstream{std::move(in)}}
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
    EXPECT_SAME_TYPE(typename it_t::value_type, char);
    EXPECT_SAME_TYPE(typename it_t::reference, char &);
    EXPECT_SAME_TYPE(typename it_t::const_reference, char &);
    EXPECT_SAME_TYPE(typename it_t::difference_type, std::ptrdiff_t);
    EXPECT_SAME_TYPE(typename it_t::size_type, size_t);
    EXPECT_SAME_TYPE(typename it_t::iterator_category, std::input_iterator_tag);
}

TEST(in_file_iterator, operations)
{
    using it_t = seqan3::detail::in_file_iterator<fake_file_t>;

    fake_file_t f{std::istringstream{"helloworld"}};

    // construct
    it_t it = f.begin();

    // deref
    EXPECT_EQ(*it, 'h');

    // pre-inc
    EXPECT_EQ(*(++it), 'e');

    // post-inc
    it++;

    // deref
    EXPECT_EQ(*it, 'l');
}

TEST(in_file_iterator, comparison)
{
    using it_t = seqan3::detail::in_file_iterator<fake_file_t>;

    fake_file_t f{std::istringstream{"helloworld"}};
    it_t it = f.begin();

    // not at end
    EXPECT_FALSE(it == std::default_sentinel);

    // consume the entire range
    ++it;
    ++it;
    ++it;
    ++it;
    ++it;
    ++it;
    ++it;
    ++it;
    ++it;
    ++it;
    ++it;

    // at end
    EXPECT_TRUE(it == std::default_sentinel);
}

TEST(in_file_iterator, file_position)
{
    using it_t = seqan3::detail::in_file_iterator<fake_file_t>;

    fake_file_t f{std::istringstream{"helloworld"}};
    it_t it = f.begin();
    auto beginning = it.file_position();

    // Go to the 6th character (w) and store it.
    ++it;
    ++it;
    ++it;
    ++it;
    ++it;
    EXPECT_EQ(*it, 'w');
    auto w_position = it.file_position();

    // Go back to the beginning.
    it.seek_to(beginning);
    EXPECT_EQ(*it, 'h');

    // Go directly to the w.
    it.seek_to(w_position);
    EXPECT_EQ(*it, 'w');
}
