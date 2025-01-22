// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iterator>

#include <seqan3/io/stream/detail/fast_istreambuf_iterator.hpp>
#include <seqan3/test/streambuf.hpp>

TEST(fast_istreambuf_iterator, concept)
{
    EXPECT_TRUE(std::input_iterator<seqan3::detail::fast_istreambuf_iterator<char>>);
    EXPECT_FALSE(std::forward_iterator<seqan3::detail::fast_istreambuf_iterator<char>>);
    EXPECT_FALSE(std::bidirectional_iterator<seqan3::detail::fast_istreambuf_iterator<char>>);
    EXPECT_FALSE(std::random_access_iterator<seqan3::detail::fast_istreambuf_iterator<char>>);
    EXPECT_FALSE((std::output_iterator<seqan3::detail::fast_istreambuf_iterator<char>, char>));
}

TEST(fast_istreambuf_iterator, construction)
{
    using type = seqan3::detail::fast_istreambuf_iterator<char>;
    EXPECT_TRUE(std::is_nothrow_default_constructible_v<type>);
    EXPECT_TRUE(std::is_nothrow_copy_constructible_v<type>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<type>);
    EXPECT_TRUE(std::is_nothrow_copy_assignable_v<type>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<type>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<type>);
    EXPECT_TRUE((std::is_constructible_v<type, std::basic_streambuf<char> &>));
}

TEST(fast_istreambuf_iterator, basic)
{
    std::istringstream str{"test"};
    seqan3::detail::fast_istreambuf_iterator<char> it{*str.rdbuf()};

    EXPECT_EQ(*it, 't');
    ++it;
    EXPECT_EQ(*it, 'e');
    it++;
    EXPECT_EQ(*it, 's');
}

TEST(fast_istreambuf_iterator, comparison)
{
    std::istringstream str{"test\n"};
    seqan3::detail::fast_istreambuf_iterator<char> it{*str.rdbuf()};

    EXPECT_FALSE(it == std::default_sentinel);
    EXPECT_FALSE(std::default_sentinel == it);
    EXPECT_TRUE(it != std::default_sentinel);
    EXPECT_TRUE(std::default_sentinel != it);
}

TEST(fast_istreambuf_iterator, cache_record_into)
{
    std::istringstream str{"record\tAAA\tBBB\tCCC\nrecord2\tXXX\tYYY\tZZZ\n"};
    seqan3::detail::fast_istreambuf_iterator<char> it{*str.rdbuf()};

    std::array<std::string_view, 4> raw_record;

    it.cache_record_into('\n', '\t', raw_record);
    ++it; // skip newline

    EXPECT_EQ(raw_record[0], "record");
    EXPECT_EQ(raw_record[1], "AAA");
    EXPECT_EQ(raw_record[2], "BBB");
    EXPECT_EQ(raw_record[3], "CCC");

    it.cache_record_into('\n', '\t', raw_record);
    ++it; // skip newline

    EXPECT_EQ(raw_record[0], "record2");
    EXPECT_EQ(raw_record[1], "XXX");
    EXPECT_EQ(raw_record[2], "YYY");
    EXPECT_EQ(raw_record[3], "ZZZ");

    EXPECT_TRUE(std::default_sentinel == it); // reached end
}

TEST(fast_istreambuf_iterator, cache_record_into_small_streambuffer)
{
    std::istringstream str{"record\tAAA\tBBB\tCCC\nrecord2\tXXX\tYYY\tZZZ\n"};
    std::istream & in{str};
    std::streambuf * orig = in.rdbuf();
    seqan3::test::streambuf_with_custom_buffer_size<3> buf(orig);
    in.rdbuf(&buf);

    seqan3::detail::fast_istreambuf_iterator<char> it{*in.rdbuf()};

    std::array<std::string_view, 4> raw_record;

    it.cache_record_into('\n', '\t', raw_record);
    ++it; // skip newline

    EXPECT_EQ(raw_record[0], "record");
    EXPECT_EQ(raw_record[1], "AAA");
    EXPECT_EQ(raw_record[2], "BBB");
    EXPECT_EQ(raw_record[3], "CCC");

    it.cache_record_into('\n', '\t', raw_record);
    ++it; // skip newline

    EXPECT_EQ(raw_record[0], "record2");
    EXPECT_EQ(raw_record[1], "XXX");
    EXPECT_EQ(raw_record[2], "YYY");
    EXPECT_EQ(raw_record[3], "ZZZ");

    EXPECT_TRUE(std::default_sentinel == it); // reached end
}

#ifndef NDEBUG
TEST(debug_fast_istreambuf_iterator, no_record_end_sign_found_after_last_field)
{
    std::istringstream str{"record\tAAA\tBBB\tCCC___oh_oh_here_is_no_newline"};
    seqan3::detail::fast_istreambuf_iterator<char> it{*str.rdbuf()};

    std::array<std::string_view, 4> raw_record;

    EXPECT_DEATH(it.cache_record_into('\n', '\t', raw_record), "");
}

TEST(debug_fast_istreambuf_iterator, not_enough_field_separation_signs_found)
{
    std::istringstream str{"record\tAAA\tBBB___oh_oh_here_is_a_tab_missing_here__CCC\n"};
    seqan3::detail::fast_istreambuf_iterator<char> it{*str.rdbuf()};

    std::array<std::string_view, 4> raw_record;

    EXPECT_DEATH(it.cache_record_into('\n', '\t', raw_record), "");
}
#endif

TEST(fast_istreambuf_iterator, cache_bytes)
{
    std::istringstream str{"ABCDEFGHIJKLMNOPQRSTUVWXYZ"};
    seqan3::detail::fast_istreambuf_iterator<char> it{*str.rdbuf()};

    std::string_view bytes = it.cache_bytes(5);
    EXPECT_EQ(bytes, "ABCDE");

    bytes = it.cache_bytes(5);
    EXPECT_EQ(bytes, "FGHIJ");
}

TEST(fast_istreambuf_iterator, cache_bytes_small_streambuffer)
{
    std::istringstream str{"ABCDEFGHIJKLMNOPQRSTUVWXYZ"};
    std::istream & in{str};
    std::streambuf * orig = in.rdbuf();
    seqan3::test::streambuf_with_custom_buffer_size<3> buf(orig);
    in.rdbuf(&buf);

    seqan3::detail::fast_istreambuf_iterator<char> it{*in.rdbuf()};

    std::string_view bytes = it.cache_bytes(5);
    EXPECT_EQ(bytes, "ABCDE");

    bytes = it.cache_bytes(5);
    EXPECT_EQ(bytes, "FGHIJ");
}

#ifndef NDEBUG
TEST(debug_fast_istreambuf_iterator, cache_bytes_too_many_bytes)
{
    std::istringstream str{"ABCDE"};
    seqan3::detail::fast_istreambuf_iterator<char> it{*str.rdbuf()};

    EXPECT_DEATH(it.cache_bytes(10), "");
}

TEST(debug_fast_istreambuf_iterator, cache_bytes_too_many_bytes_small_streambuffer)
{
    std::istringstream str{"ABCDE"};
    std::istream & in{str};
    std::streambuf * orig = in.rdbuf();
    seqan3::test::streambuf_with_custom_buffer_size<3> buf(orig);
    in.rdbuf(&buf);

    seqan3::detail::fast_istreambuf_iterator<char> it{*in.rdbuf()};

    EXPECT_DEATH(it.cache_bytes(10), "");
}
#endif
