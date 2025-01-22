// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <forward_list>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/utility/views/repeat.hpp>
#include <seqan3/utility/views/zip.hpp>

#include "../../range/iterator_test_template.hpp"

using range_t = std::vector<int>;
using other_range_t = std::vector<std::string>;
using zip_view_t = decltype(seqan3::views::zip(std::declval<range_t &>(), std::declval<other_range_t &>()));
using zip_iterator_t = std::ranges::iterator_t<zip_view_t>;

template <>
struct iterator_fixture<zip_iterator_t> : public ::testing::Test
{
    using iterator_tag = std::random_access_iterator_tag;

    static constexpr bool const_iterable = true;

    range_t range{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    other_range_t other_range{"AA", "BBB", "CC", "DDD"};

    std::vector<seqan3::common_pair<int, std::string>> expected_range{{0, "AA"}, {1, "BBB"}, {2, "CC"}, {3, "DDD"}};

    zip_view_t test_range{seqan3::views::zip(range, other_range)};
};
INSTANTIATE_TYPED_TEST_SUITE_P(zip_iterator_test, iterator_fixture, zip_iterator_t, );

using seqan3::operator""_dna4;

class zip_test : public ::testing::Test
{
protected:
    using range_t = std::vector<int>;
    using const_range_t = std::vector<int> const;
    using other_range_t = std::vector<std::string>;
    using forward_range_t = std::forward_list<int>;
    using view_t = decltype(seqan3::views::repeat('L'));

    range_t range{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    const_range_t const_range{range};
    other_range_t other_range{"AA", "BBB", "CC", "DDD"};
    forward_range_t forward_range{range.begin(), range.end()};
    static constexpr view_t view{};

    using zip_common_range_t = decltype(seqan3::views::zip(range, other_range));
    using zip_not_common_range_t = decltype(seqan3::views::zip(range, other_range, view));
    using zip_const_range_t = decltype(seqan3::views::zip(range, const_range));
    using zip_forward_range_t = decltype(seqan3::views::zip(range, other_range, forward_range));
    using const_zip_t = decltype(seqan3::views::zip(range, other_range)) const;

    zip_common_range_t zip_common_range()
    {
        return seqan3::views::zip(range, other_range);
    }

    zip_not_common_range_t zip_not_common_range()
    {
        return seqan3::views::zip(range, other_range, view);
    }

    zip_const_range_t zip_const_range()
    {
        return seqan3::views::zip(range, const_range);
    }

    zip_forward_range_t zip_forward_range()
    {
        return seqan3::views::zip(range, other_range, forward_range);
    }

    const_zip_t const_zip()
    {
        return seqan3::views::zip(range, other_range);
    }
};

TEST_F(zip_test, concepts)
{
    EXPECT_TRUE(std::ranges::forward_range<zip_forward_range_t>);
    EXPECT_FALSE(std::ranges::bidirectional_range<zip_forward_range_t>);

    EXPECT_TRUE(std::ranges::random_access_range<zip_common_range_t>);
    EXPECT_TRUE(std::ranges::random_access_range<zip_not_common_range_t>);
    EXPECT_TRUE(std::ranges::random_access_range<zip_const_range_t>);
    EXPECT_TRUE(std::ranges::random_access_range<const_zip_t>);

    EXPECT_TRUE(std::ranges::view<zip_common_range_t>);
    EXPECT_TRUE(std::ranges::view<zip_not_common_range_t>);
    EXPECT_TRUE(std::ranges::view<zip_const_range_t>);
    EXPECT_TRUE(std::ranges::view<zip_forward_range_t>);
    EXPECT_FALSE(std::ranges::view<const_zip_t>); // const lvalue is not movable, because it cannot be assigned to

    EXPECT_TRUE(std::ranges::sized_range<zip_common_range_t>);
    EXPECT_FALSE(std::ranges::sized_range<zip_not_common_range_t>); // seqan3::views::repeat has no size (infinite)
    EXPECT_TRUE(std::ranges::sized_range<zip_const_range_t>);
    EXPECT_FALSE(std::ranges::sized_range<zip_forward_range_t>); // std::forward_list is not sized
    EXPECT_TRUE(std::ranges::sized_range<const_zip_t>);

    EXPECT_TRUE(std::ranges::common_range<zip_common_range_t>);
    EXPECT_FALSE(std::ranges::common_range<zip_not_common_range_t>);
    EXPECT_TRUE(std::ranges::common_range<zip_const_range_t>);
    EXPECT_TRUE(std::ranges::common_range<zip_forward_range_t>);
    EXPECT_TRUE(std::ranges::common_range<const_zip_t>);

    EXPECT_TRUE(std::ranges::common_range<zip_common_range_t>);
    EXPECT_FALSE(std::ranges::common_range<zip_not_common_range_t>);
    EXPECT_TRUE(std::ranges::common_range<zip_const_range_t>);
    EXPECT_TRUE(std::ranges::common_range<zip_forward_range_t>);
    EXPECT_TRUE(std::ranges::common_range<const_zip_t>);

    EXPECT_TRUE((std::ranges::output_range<zip_common_range_t, std::pair<int &, std::string &>>));
    EXPECT_TRUE((std::ranges::output_range<zip_not_common_range_t, std::tuple<int &, std::string &, char &>>));
    EXPECT_FALSE((std::ranges::output_range<zip_const_range_t, std::pair<int &, int &>>));
    EXPECT_FALSE((std::ranges::output_range<zip_const_range_t, std::pair<int &, int const &>>));
    EXPECT_TRUE((std::ranges::output_range<zip_forward_range_t, std::tuple<int &, std::string &, int &>>));
    EXPECT_TRUE((std::ranges::output_range<const_zip_t, std::pair<int &, std::string &>>));
}

TEST_F(zip_test, basic)
{
    {
        auto zip_view = zip_common_range();
        size_t i{};
        for (auto && [elem_1, elem_2] : zip_view)
        {
            EXPECT_EQ(elem_1, range[i]);
            EXPECT_EQ(elem_2, other_range[i]);
            ++i;
        }
        EXPECT_EQ(i, 4u);
        EXPECT_EQ(zip_view.size(), 4u);
    }
    {
        auto zip_view = seqan3::views::zip(range, other_range);
        size_t i{};
        for (auto && [elem_1, elem_2] : zip_view)
        {
            EXPECT_EQ(elem_1, range[i]);
            EXPECT_EQ(elem_2, other_range[i]);
            ++i;
        }
        EXPECT_EQ(i, 4u);
        EXPECT_EQ(zip_view.size(), 4u);
    }
}

TEST_F(zip_test, combine)
{
    auto zip_view = zip_common_range() | std::views::take(2);
    size_t i{};
    for (auto && [elem_1, elem_2] : zip_view)
    {
        EXPECT_EQ(elem_1, range[i]);
        EXPECT_EQ(elem_2, other_range[i]);
        ++i;
    }
    EXPECT_EQ(i, 2u);
    EXPECT_EQ(zip_view.size(), 2u);
}

TEST_F(zip_test, alignment_usage_1)
{
    seqan3::dna4_vector sequence_1{"AAAAA"_dna4};
    seqan3::dna4_vector sequence_2{"TTTTT"_dna4};
    std::tuple<seqan3::dna4_vector, seqan3::dna4_vector> sequence_pair{sequence_1, sequence_2};
    auto tuple_view = std::views::single(sequence_pair);
    auto zipped_tuple = seqan3::views::zip(tuple_view, std::views::iota(0));
    auto chunked_zip = zipped_tuple | seqan3::views::chunk(1);
    (void)chunked_zip;
}

TEST_F(zip_test, alignment_usage_2)
{
    seqan3::dna4_vector sequence_1{"AAAAA"_dna4};
    seqan3::dna4_vector sequence_2{"TTTTT"_dna4};
    auto tuple_view = std::views::single(std::tie(sequence_1, sequence_2));
    auto zipped_tuple = seqan3::views::zip(tuple_view, std::views::iota(0));
    auto chunked_zip = zipped_tuple | seqan3::views::chunk(1);
    (void)chunked_zip;
}

TEST_F(zip_test, use_as_output_range)
{
    auto zip_view = zip_common_range();
    *zip_view.begin() = std::pair(23, "FF");
    EXPECT_EQ(std::get<0>(*zip_view.begin()), 23);
    EXPECT_EQ(std::get<1>(*zip_view.begin()), "FF");

    size_t i{1u};
    for (auto && [elem_1, elem_2] : zip_view | std::views::drop(1))
    {
        EXPECT_EQ(elem_1, range[i]);
        EXPECT_EQ(elem_2, other_range[i]);
        ++i;
    }
    EXPECT_EQ(i, 4u);
    EXPECT_EQ(zip_view.size(), 4u);
}
