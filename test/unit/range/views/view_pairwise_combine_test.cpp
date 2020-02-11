// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <forward_list>
#include <list>
#include <utility>
#include <vector>

#include <range/v3/algorithm/equal.hpp>
#include <range/v3/view/filter.hpp>

#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/views/pairwise_combine.hpp>
#include <seqan3/range/views/take.hpp>
#include <seqan3/test/pretty_printing.hpp>

template <typename t>
class pairwise_combine_base_test : public ::testing::Test
{
public:

    using view_t       = decltype(seqan3::detail::pairwise_combine_view{std::views::all(std::declval<t &>())});
    using const_view_t = decltype(seqan3::detail::pairwise_combine_view{
                                    std::views::all(std::declval<t const &>())});

    auto create_view()
    {
        return view_t{container};
    }

    t & resource()
    {
        return container;
    }

    auto const & expect()
    {
        return res;
    }

protected:

    void SetUp() override
    {
        if constexpr (seqan3::detail::is_type_specialisation_of_v<t, std::forward_list>)
        {
            container.push_front('d');
            container.push_front('c');
            container.push_front('b');
            container.push_front('a');
        }
        else
        {
            container.push_back('a');
            container.push_back('b');
            container.push_back('c');
            container.push_back('d');
        }
        res.push_back(std::tuple{'a', 'b'});
        res.push_back(std::tuple{'a', 'c'});
        res.push_back(std::tuple{'a', 'd'});
        res.push_back(std::tuple{'b', 'c'});
        res.push_back(std::tuple{'b', 'd'});
        res.push_back(std::tuple{'c', 'd'});
    }

private:
    t container;

    std::vector<std::tuple<char, char>> res;
};

template <typename t>
class pairwise_combine_test : public pairwise_combine_base_test<t>
{};

template <typename t>
class pairwise_combine_iterator_test : public pairwise_combine_base_test<t>
{};

using test_types = ::testing::Types<std::vector<char>, std::list<char>, std::forward_list<char>>;

TYPED_TEST_SUITE(pairwise_combine_test, test_types, );
TYPED_TEST_SUITE(pairwise_combine_iterator_test, test_types, );

TYPED_TEST(pairwise_combine_iterator_test, concepts)
{
    EXPECT_TRUE(std::forward_iterator<std::ranges::iterator_t<typename TestFixture::view_t>>);

    if constexpr (std::bidirectional_iterator<std::ranges::iterator_t<TypeParam>>)
    {
        EXPECT_TRUE(std::bidirectional_iterator<std::ranges::iterator_t<typename TestFixture::view_t>>);
    }

    if constexpr (std::random_access_iterator<std::ranges::iterator_t<TypeParam>>)
    {
        EXPECT_TRUE(std::random_access_iterator<std::ranges::iterator_t<typename TestFixture::view_t>>);
    }
}

TYPED_TEST(pairwise_combine_iterator_test, construction)
{
    using iterator_t = std::ranges::iterator_t<typename TestFixture::view_t>;
    using underlying_iter_t = std::ranges::iterator_t<TypeParam>;

    EXPECT_TRUE(std::is_nothrow_default_constructible_v<iterator_t>);
    EXPECT_TRUE(std::is_nothrow_copy_constructible_v<iterator_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<iterator_t>);
    EXPECT_TRUE(std::is_nothrow_copy_assignable_v<iterator_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<iterator_t>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<iterator_t>);
    EXPECT_TRUE((std::is_nothrow_constructible_v<iterator_t, underlying_iter_t, underlying_iter_t, underlying_iter_t>));

    using const_itertor_t = std::ranges::iterator_t<typename TestFixture::view_t const>;
    EXPECT_TRUE((std::is_nothrow_constructible_v<const_itertor_t, iterator_t>));

    auto r = this->resource();
    iterator_t it{r.begin(), r.begin(), r.end()};

    const_itertor_t cit{it};
    EXPECT_TRUE(it == cit);

    const_itertor_t cit2{};
    cit2 = it;
    EXPECT_TRUE(cit == cit2);
}

TYPED_TEST(pairwise_combine_iterator_test, associated_types)
{
    // The iterator over the view over the non-const range
    using v_iter_t        = std::ranges::iterator_t<typename TestFixture::view_t>;
    // The iterator over the view over the const range
    using v_const_iter_t  = std::ranges::iterator_t<typename TestFixture::const_view_t>;
    // The non-const range iterator
    using u_iter_t        = std::ranges::iterator_t<TypeParam>;
    // The const range iterator
    using u_const_iter_t  = std::ranges::iterator_t<TypeParam const>;

    { // non-const view over non-const range
        using u_ref_t = typename std::iterator_traits<u_iter_t>::reference;
        EXPECT_TRUE((std::is_same_v<typename std::iterator_traits<v_iter_t>::reference,
                                    seqan3::common_tuple<u_ref_t, u_ref_t>>));
    }

    { // non-const view over const range
        using u_ref_t = typename std::iterator_traits<u_const_iter_t>::reference;
        EXPECT_TRUE((std::is_same_v<typename std::iterator_traits<v_const_iter_t>::reference,
                                    seqan3::common_tuple<u_ref_t, u_ref_t>>));
    }

    { // const view over non-const range
        using u_ref_t = typename std::iterator_traits<u_iter_t>::reference;
        EXPECT_TRUE((std::is_same_v<typename std::iterator_traits<const v_iter_t>::reference,
                                    seqan3::common_tuple<u_ref_t, u_ref_t>>));
    }

    { // const view over const range
        using u_ref_t = typename std::iterator_traits<u_const_iter_t>::reference;
        EXPECT_TRUE((std::is_same_v<typename std::iterator_traits<const v_const_iter_t>::reference,
                                    seqan3::common_tuple<u_ref_t, u_ref_t>>));
    }

    { // value type
        using u_val_t = typename std::iterator_traits<u_iter_t>::value_type;
        EXPECT_TRUE((std::is_same_v<typename std::iterator_traits<v_iter_t>::value_type,
                                    std::tuple<u_val_t, u_val_t>>));
    }
}

TYPED_TEST(pairwise_combine_iterator_test, dereference)
{
    using iterator_t = std::ranges::iterator_t<typename TestFixture::view_t>;

    auto r = this->resource();
    iterator_t it{r.begin(), r.begin(), r.end()};

    EXPECT_EQ(*it, (std::tuple{'a', 'b'}));

    std::get<0>(*it) = 'x';

    iterator_t const c_it{it};
    EXPECT_EQ(*c_it, (std::tuple{'x', 'b'}));
}

TYPED_TEST(pairwise_combine_iterator_test, pre_increment)
{
    using iterator_t = std::ranges::iterator_t<typename TestFixture::view_t>;

    auto r = this->resource();
    iterator_t it{r.begin(), r.begin(), r.end()};

    EXPECT_EQ(*(++it), (std::tuple{'a', 'c'}));
}

TYPED_TEST(pairwise_combine_iterator_test, post_increment)
{
    using iterator_t = std::ranges::iterator_t<typename TestFixture::view_t>;

    auto r = this->resource();
    iterator_t it{r.begin(), r.begin(), r.end()};

    EXPECT_EQ(*(it++), (std::tuple{'a', 'b'}));
    EXPECT_EQ(*it, (std::tuple{'a', 'c'}));
}

TYPED_TEST(pairwise_combine_iterator_test, pre_decrement)
{
    if constexpr (std::ranges::bidirectional_range<TypeParam>)
    {
        using iterator_t = std::ranges::iterator_t<typename TestFixture::view_t>;

        auto r = this->resource();
        iterator_t it{std::ranges::prev(r.end()), r.begin(), r.end()};

        if constexpr (std::bidirectional_iterator<decltype(r.begin())>)
        {
            EXPECT_EQ(*(--it), (std::tuple{'c', 'd'}));
        }
    }
}

TYPED_TEST(pairwise_combine_iterator_test, post_decrement)
{
    if constexpr (std::ranges::bidirectional_range<TypeParam>)
    {
        using iterator_t = std::ranges::iterator_t<typename TestFixture::view_t>;

        auto r = this->resource();
        iterator_t it{std::ranges::prev(r.end()), r.begin(), r.end()};

        if constexpr (std::bidirectional_iterator<decltype(r.begin())>)
        {
            --it;
            EXPECT_EQ(*(it--), (std::tuple{'c', 'd'}));
            EXPECT_EQ(*it, (std::tuple{'b', 'd'}));
        }
    }
}

TYPED_TEST(pairwise_combine_iterator_test, equality)
{
    using iterator_t = std::ranges::iterator_t<typename TestFixture::view_t>;

    auto r = this->resource();
    iterator_t it{r.begin(), r.begin(), r.end()};

    iterator_t it_2 = it++;
    EXPECT_TRUE(it == it);
    EXPECT_TRUE(it != it_2);
}

TYPED_TEST(pairwise_combine_iterator_test, subscript)
{
    using iterator_t = std::ranges::iterator_t<typename TestFixture::view_t>;

    auto r = this->resource();
    iterator_t it{r.begin(), r.begin(), r.end()};

    if constexpr (std::random_access_iterator<decltype(r.begin())>)
    {
        EXPECT_EQ(it[0], (std::tuple{'a', 'b'}));
        EXPECT_EQ(it[1], (std::tuple{'a', 'c'}));
        EXPECT_EQ(it[2], (std::tuple{'a', 'd'}));
        EXPECT_EQ(it[3], (std::tuple{'b', 'c'}));
        EXPECT_EQ(it[4], (std::tuple{'b', 'd'}));
        EXPECT_EQ(it[5], (std::tuple{'c', 'd'}));
    }
}

TYPED_TEST(pairwise_combine_iterator_test, advance_n)
{
    using iterator_t = std::ranges::iterator_t<typename TestFixture::view_t>;

    auto r = this->resource();
    iterator_t it{r.begin(), r.begin(), r.end()};

    if constexpr (std::random_access_iterator<decltype(r.begin())>)
    {
        it = it + 1;
        EXPECT_EQ(*it, (std::tuple{'a', 'c'}));
        it += 2;
        EXPECT_EQ(*it, (std::tuple{'b', 'c'}));
        it = 2 + it;
        EXPECT_EQ(*it, (std::tuple{'c', 'd'}));
    }
}

TYPED_TEST(pairwise_combine_iterator_test, decrement_n)
{
    if constexpr (std::ranges::random_access_range<TypeParam>)
    {
        using iterator_t = std::ranges::iterator_t<typename TestFixture::view_t>;

        auto r = this->resource();
        iterator_t it{std::ranges::prev(r.end()), r.begin(), r.end()};

        it = it - 2;
        EXPECT_EQ(*it, (std::tuple{'b', 'd'}));
        it -= 3;
        EXPECT_EQ(*it, (std::tuple{'a', 'c'}));
    }
}

TYPED_TEST(pairwise_combine_iterator_test, distance)
{
    if constexpr (std::ranges::random_access_range<TypeParam>)
    {
        using iterator_t = std::ranges::iterator_t<typename TestFixture::view_t>;

        auto r = this->resource();
        iterator_t it_1{r.begin(), r.begin(), r.end()};
        iterator_t it_2{std::ranges::prev(r.end()), r.begin(), r.end()};

        EXPECT_EQ(it_2 - it_1, 6);
        --it_2;
        EXPECT_EQ(it_2 - it_1, 5);
        ++it_1;
        EXPECT_EQ(it_2 - it_1, 4);
    }
}

TYPED_TEST(pairwise_combine_iterator_test, order)
{
    using iterator_t = std::ranges::iterator_t<typename TestFixture::view_t>;

    auto r = this->resource();
    iterator_t it_1{r.begin(), r.begin(), r.end()};
    iterator_t it_2{r.begin(), r.begin(), r.end()};
    ++it_2;

    if constexpr (std::totally_ordered<decltype(r.begin())>)
    {
        EXPECT_FALSE(it_1 < it_1);
        EXPECT_TRUE(it_1 <= it_1);
        EXPECT_FALSE(it_1 > it_1);
        EXPECT_TRUE(it_1 >= it_1);

        EXPECT_TRUE(it_1 < it_2);
        EXPECT_TRUE(it_1 <= it_2);
        EXPECT_FALSE(it_1 > it_2);
        EXPECT_FALSE(it_1 >= it_2);
    }
}

TYPED_TEST(pairwise_combine_test, view_concept)
{
    EXPECT_TRUE(std::ranges::input_range<typename TestFixture::view_t>);
    EXPECT_TRUE(std::ranges::forward_range<typename TestFixture::view_t>);
    EXPECT_TRUE(std::ranges::view<typename TestFixture::view_t>);
    EXPECT_TRUE((std::ranges::output_range<typename TestFixture::view_t, std::tuple<char &, char &>>));
    EXPECT_EQ(std::ranges::bidirectional_range<TypeParam>,
              std::ranges::bidirectional_range<typename TestFixture::view_t>);
    EXPECT_EQ(std::ranges::sized_range<TypeParam>, std::ranges::sized_range<typename TestFixture::view_t>);
    EXPECT_EQ(std::ranges::random_access_range<TypeParam>, std::ranges::random_access_range<typename TestFixture::view_t>);

    EXPECT_TRUE(std::ranges::input_range<typename TestFixture::const_view_t>);
    EXPECT_TRUE(std::ranges::forward_range<typename TestFixture::const_view_t>);
    EXPECT_TRUE(std::ranges::view<typename TestFixture::const_view_t>);
    //TODO: Investigate, because this should be false. It cannot be used as a output range!
    // EXPECT_FALSE((std::ranges::output_range<typename TestFixture::const_view_t, std::tuple<char, char>>));
    EXPECT_EQ(std::ranges::bidirectional_range<TypeParam>,
              std::ranges::bidirectional_range<typename TestFixture::const_view_t>);
    EXPECT_EQ(std::ranges::sized_range<TypeParam>, std::ranges::sized_range<typename TestFixture::const_view_t>);
    EXPECT_EQ(std::ranges::random_access_range<TypeParam>, std::ranges::random_access_range<typename TestFixture::const_view_t>);
}

TYPED_TEST(pairwise_combine_test, basic_construction)
{
    using view_t = typename TestFixture::view_t;

    EXPECT_TRUE(std::is_default_constructible_v<view_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<view_t>);
    EXPECT_TRUE(std::is_move_constructible_v<view_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<view_t>);
    EXPECT_TRUE(std::is_move_assignable_v<view_t>);
    EXPECT_TRUE(std::is_destructible_v<view_t>);
}

TYPED_TEST(pairwise_combine_test, begin)
{
    auto v = this->create_view();
    auto const cv{v};

    EXPECT_EQ(*v.begin(), (std::tuple{'a', 'b'}));
    EXPECT_EQ(*cv.begin(), (std::tuple{'a', 'b'}));
    EXPECT_EQ(*v.cbegin(), (std::tuple{'a', 'b'}));
}

TYPED_TEST(pairwise_combine_test, end)
{
    auto v = this->create_view();
    auto const cv{v};

    EXPECT_TRUE(v.begin() != v.end());
    EXPECT_TRUE(cv.begin() != cv.end());
    EXPECT_TRUE(v.cbegin() != v.cend());
}

TYPED_TEST(pairwise_combine_test, iterate)
{
    auto v = this->create_view();
    using ref_t = std::iter_reference_t<std::ranges::iterator_t<decltype(v)>>;
    std::vector<ref_t> cmp;

    for (auto r : v)
        cmp.push_back(r);

    EXPECT_TRUE(ranges::equal(cmp, this->expect()));
}

TYPED_TEST(pairwise_combine_test, iterate_reverse)
{
    if constexpr (std::ranges::bidirectional_range<TypeParam>)
    {
        auto v = this->create_view();
        using ref_t = std::iter_reference_t<std::ranges::iterator_t<decltype(v)>>;
        std::vector<ref_t> cmp;

        for (auto r : v | std::views::reverse)
            cmp.push_back(r);

        EXPECT_TRUE(ranges::equal(cmp, this->expect() | std::views::reverse));
    }
}

TYPED_TEST(pairwise_combine_test, size)
{
    auto v = this->create_view();

    if constexpr (std::ranges::sized_range<decltype(v)>)
    {
        EXPECT_EQ(std::ranges::size(v), 6u);
        EXPECT_EQ(v.size(), 6u);
    }
}

TEST(pairwise_combine_fn_test, filter_output)
{
    std::vector orig{'a', 'b', 'x', 'c', 'd'};

    auto v = orig | seqan3::views::pairwise_combine;

    using ref_t = std::iter_reference_t<std::ranges::iterator_t<decltype(v)>>;
    std::vector<ref_t> cmp;

    auto v_filter = v | std::views::filter([](auto tpl)
                        {
                            return !((std::get<0>(tpl) == 'x') || (std::get<1>(tpl) == 'x'));
                        });

    for (auto r : v_filter)
        cmp.push_back(r);

    auto it = std::ranges::begin(cmp);
    EXPECT_EQ(*it, (std::tuple{'a', 'b'}));
    EXPECT_EQ(*++it, (std::tuple{'a', 'c'}));
    EXPECT_EQ(*++it, (std::tuple{'a', 'd'}));
    EXPECT_EQ(*++it, (std::tuple{'b', 'c'}));
    EXPECT_EQ(*++it, (std::tuple{'b', 'd'}));
    EXPECT_EQ(*++it, (std::tuple{'c', 'd'}));
}

TEST(pairwise_combine_fn_test, filter_input)
{
    std::vector orig{'a', 'b', 'x', 'c', 'd'};
    auto v_filter = orig | std::views::filter([](char c) { return c != 'x'; });

    auto v = v_filter | seqan3::views::pairwise_combine;

    using ref_t = std::iter_reference_t<std::ranges::iterator_t<decltype(v)>>;
    std::vector<ref_t> cmp;

    for (auto r : v)
        cmp.push_back(r);

    auto it = std::ranges::begin(cmp);
    EXPECT_EQ(*it, (std::tuple{'a', 'b'}));
    EXPECT_EQ(*++it, (std::tuple{'a', 'c'}));
    EXPECT_EQ(*++it, (std::tuple{'a', 'd'}));
    EXPECT_EQ(*++it, (std::tuple{'b', 'c'}));
    EXPECT_EQ(*++it, (std::tuple{'b', 'd'}));
    EXPECT_EQ(*++it, (std::tuple{'c', 'd'}));
}

TEST(pairwise_combine_fn_test, output)
{
    std::vector orig{'a', 'b', 'c', 'd'};
    auto v = orig | seqan3::views::pairwise_combine;

    *v.begin() = std::tuple{'x', 'y'};

    auto it = v.begin();
    EXPECT_EQ(*it,   (std::tuple{'x', 'y'}));
    EXPECT_EQ(*++it, (std::tuple{'x', 'c'}));
    EXPECT_EQ(*++it, (std::tuple{'x', 'd'}));
    EXPECT_EQ(*++it, (std::tuple{'y', 'c'}));
    EXPECT_EQ(*++it, (std::tuple{'y', 'd'}));
    EXPECT_EQ(*++it, (std::tuple{'c', 'd'}));
}

TEST(pairwise_combine_fn_test, const_source)
{
    std::vector orig{'a', 'b', 'c', 'd'};
    std::vector<char> const c_orig{orig};
    auto v = c_orig | seqan3::views::pairwise_combine;

    using ref_t = typename std::iterator_traits<std::ranges::iterator_t<decltype(v)>>::reference;
    std::vector<ref_t> cmp;

    for (auto r : v)
        cmp.push_back(r);

    auto it = std::ranges::begin(cmp);
    EXPECT_EQ(*it, (std::tuple{'a', 'b'}));
    EXPECT_EQ(*++it, (std::tuple{'a', 'c'}));
    EXPECT_EQ(*++it, (std::tuple{'a', 'd'}));
    EXPECT_EQ(*++it, (std::tuple{'b', 'c'}));
    EXPECT_EQ(*++it, (std::tuple{'b', 'd'}));
    EXPECT_EQ(*++it, (std::tuple{'c', 'd'}));
}
