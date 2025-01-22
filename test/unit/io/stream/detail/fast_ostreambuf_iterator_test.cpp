// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <iterator>

#include <seqan3/io/stream/detail/fast_ostreambuf_iterator.hpp>
#include <seqan3/test/expect_same_type.hpp>

TEST(fast_ostreambuf_iterator, concept)
{
    using type = seqan3::detail::fast_ostreambuf_iterator<char>;

    EXPECT_TRUE((std::output_iterator<type, char>));

    EXPECT_FALSE(std::input_iterator<type>);
    EXPECT_FALSE(std::forward_iterator<type>);
    EXPECT_FALSE(std::bidirectional_iterator<type>);
    EXPECT_FALSE(std::random_access_iterator<type>);
}

TEST(fast_ostreambuf_iterator, construction)
{
    using type = seqan3::detail::fast_ostreambuf_iterator<char>;

    EXPECT_TRUE(std::is_nothrow_default_constructible_v<type>);
    EXPECT_TRUE(std::is_nothrow_copy_constructible_v<type>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<type>);
    EXPECT_TRUE(std::is_nothrow_copy_assignable_v<type>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<type>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<type>);
    EXPECT_TRUE((std::is_constructible_v<type, std::basic_streambuf<char> &>));
}

TEST(fast_ostreambuf_iterator, assignment)
{
    std::ostringstream ostr{};
    seqan3::detail::fast_ostreambuf_iterator<char> it{*ostr.rdbuf()};

    it = 't';
    *it = 'e';

    it++; // no op
    ++it; // no op

    *it++ = 's';
    *++it = 't';

    EXPECT_EQ(ostr.str(), std::string_view{"test"});
}

TEST(fast_ostreambuf_iterator, write_range_simple_case)
{
    // the simple case without using the return value of write_value
    std::ostringstream ostr{};
    seqan3::detail::fast_ostreambuf_iterator<char> it{*ostr.rdbuf()};

    std::string rng{"test\ntestest"};

    it.write_range(rng);

    EXPECT_EQ(ostr.str(), std::string_view{"test\ntestest"});
}

TEST(fast_ostreambuf_iterator, write_range_ensure_overflow)
{
    // ensure that buffer overflows at least once
    std::ostringstream ostr;
    char buf[40];
    ostr.rdbuf()->pubsetbuf(buf, sizeof buf);

    seqan3::detail::fast_ostreambuf_iterator<char> it{*ostr.rdbuf()};

    std::string rng{"veryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryverylargerange"};

    it.write_range(rng);

    EXPECT_EQ(ostr.str().substr(0, 78),
              std::string_view{"veryveryveryveryveryveryveryveryveryveryveryveryveryveryveryveryverylargerange"});
}

TEST(fast_ostreambuf_iterator, write_range_return_value)
{
    // use the return value of write_range to keep track of the chunk you have already written
    std::ostringstream ostr{};
    seqan3::detail::fast_ostreambuf_iterator<char> it{*ostr.rdbuf()};

    std::string rng{"test\ntestest"};
    using str_it_t = std::string::iterator;

    auto rng_it = it.write_range(std::ranges::subrange<str_it_t, str_it_t>(rng.begin(), rng.begin() + 5));

    EXPECT_EQ(rng_it, (rng.begin() + 5));
    EXPECT_EQ(ostr.str(), std::string_view{"test\n"});

    rng_it = it.write_range(std::ranges::subrange<str_it_t, str_it_t>(rng.begin() + 5, rng.end()));

    EXPECT_EQ(rng_it, rng.end());
    EXPECT_EQ(ostr.str(), std::string_view{"test\ntestest"});
}

TEST(fast_ostreambuf_iterator, write_range_unsafe_range)
{
    // write_range on a non-forwarding range returns void
    std::ostringstream ostr{};
    seqan3::detail::fast_ostreambuf_iterator<char> it{*ostr.rdbuf()};

    EXPECT_SAME_TYPE(void, decltype(it.write_range(std::string{"foo"})));

    it.write_range(std::string{"foo"});

    EXPECT_EQ(ostr.str(), std::string_view{"foo"});
}

TEST(fast_ostreambuf_iterator, write_number)
{
    {
        std::ostringstream ostr;
        char buf[400];
        ostr.rdbuf()->pubsetbuf(buf, sizeof buf);

        seqan3::detail::fast_ostreambuf_iterator<char> it{*ostr.rdbuf()};

        uint64_t num{54'389'234};

        it.write_number(num);

        EXPECT_EQ(ostr.str().substr(0, 8), std::string_view{"54389234"});
    }

    {
        std::ostringstream ostr;
        char buf[100];
        ostr.rdbuf()->pubsetbuf(buf, sizeof buf);

        seqan3::detail::fast_ostreambuf_iterator<char> it{*ostr.rdbuf()};

        uint64_t num{54'389'234};

        it.write_number(num);

        EXPECT_EQ(ostr.str().substr(0, 8), std::string_view{"54389234"});
    }
}

TEST(fast_ostreambuf_iterator, write_end_of_line)
{
    std::ostringstream ostr{};
    seqan3::detail::fast_ostreambuf_iterator<char> it{*ostr.rdbuf()};

    ostr << "test";
    it.write_end_of_line(false);

    ostr << "testest";
    it.write_end_of_line(true);

    EXPECT_EQ(ostr.str(), std::string_view{"test\ntestest\r\n"});
}
