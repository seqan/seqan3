// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>
#include <vector>

#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/utility/views/enforce_random_access.hpp>

#include "../../range/iterator_test_template.hpp"

class common_pseudo_random_access_range
{
public:
    using urng_t = std::vector<int>;

    common_pseudo_random_access_range() = default;
    common_pseudo_random_access_range(urng_t urng) : urng{std::move(urng)}
    {}

    template <typename u_iterator_t>
    class test_iterator : public seqan3::detail::inherited_iterator_base<test_iterator<u_iterator_t>, u_iterator_t>
    {
    public:
        using base_t = seqan3::detail::inherited_iterator_base<test_iterator<u_iterator_t>, u_iterator_t>;
        using iterator_category = std::bidirectional_iterator_tag;
        using iterator_concept = iterator_category;

        using base_t::base_t;
    };

    auto begin() noexcept
    {
        return test_iterator<typename urng_t::iterator>{urng.begin()};
    }

    auto begin() const
    {
        return test_iterator<typename urng_t::const_iterator>{urng.begin()};
    }

    auto end() noexcept
    {
        return test_iterator<typename urng_t::iterator>{urng.end()};
    }

    auto end() const
    {
        return test_iterator<typename urng_t::const_iterator>{urng.end()};
    }

    std::vector<int> urng{};
};

class sentinel_pseudo_random_access_range : public common_pseudo_random_access_range
{
public:
    template <typename u_iterator_t>
    class test_iterator : public seqan3::detail::inherited_iterator_base<test_iterator<u_iterator_t>, u_iterator_t>
    {
    public:
        using base_t = seqan3::detail::inherited_iterator_base<test_iterator<u_iterator_t>, u_iterator_t>;
        using iterator_category = std::bidirectional_iterator_tag;
        using iterator_concept = iterator_category;

        using base_t::base_t;

        test_iterator(u_iterator_t it, u_iterator_t last) : base_t{it}, last{last}
        {}

        using base_t::operator==;
        using base_t::operator!=;
        using base_t::operator-;
        bool operator==(std::default_sentinel_t const &) const
        {
            return this->base() == last;
        }

        friend bool operator==(std::default_sentinel_t const &, test_iterator const & rhs)
        {
            return rhs == std::default_sentinel;
        }

        bool operator!=(std::default_sentinel_t const &) const
        {
            return !(*this == std::default_sentinel);
        }

        friend bool operator!=(std::default_sentinel_t const &, test_iterator const & rhs)
        {
            return rhs != std::default_sentinel;
        }

        typename base_t::difference_type operator-(std::default_sentinel_t const &) const
        {
            return this->base() - this->last;
        }

        friend typename base_t::difference_type operator-(std::default_sentinel_t const &, test_iterator const & rhs)
        {
            return rhs.last - rhs.base();
        }

        u_iterator_t last{};
    };

    auto begin() noexcept
    {
        return test_iterator<typename urng_t::iterator>{urng.begin(), urng.end()};
    }

    auto begin() const
    {
        return test_iterator<typename urng_t::const_iterator>{urng.begin(), urng.end()};
    }

    auto end() noexcept
    {
        return std::default_sentinel;
    }

    auto end() const
    {
        return std::default_sentinel;
    }
};

template <typename t>
class enforce_random_access_test : public ::testing::Test
{};

using testing_types =
    ::testing::Types<std::vector<int>, common_pseudo_random_access_range, sentinel_pseudo_random_access_range>;

TYPED_TEST_SUITE(enforce_random_access_test, testing_types, );

TYPED_TEST(enforce_random_access_test, concepts)
{
    using enforce_random_access_type = decltype(std::declval<TypeParam &>() | seqan3::views::enforce_random_access);

    // guaranteed concepts
    EXPECT_TRUE(std::ranges::random_access_range<enforce_random_access_type>);
    EXPECT_TRUE(std::ranges::view<enforce_random_access_type>);
    EXPECT_TRUE(std::ranges::viewable_range<enforce_random_access_type>);

    // preserved concepts
    EXPECT_EQ(std::ranges::sized_range<TypeParam>, std::ranges::sized_range<enforce_random_access_type>);
    EXPECT_EQ(std::ranges::common_range<TypeParam>, std::ranges::common_range<enforce_random_access_type>);
    EXPECT_EQ(std::ranges::contiguous_range<TypeParam>, std::ranges::contiguous_range<enforce_random_access_type>);
    EXPECT_EQ(seqan3::const_iterable_range<TypeParam>, seqan3::const_iterable_range<enforce_random_access_type>);
    EXPECT_EQ((std::ranges::output_range<TypeParam, int>),
              (std::ranges::output_range<enforce_random_access_type, int>));
}

TYPED_TEST(enforce_random_access_test, adaptor)
{

    std::vector<int> source{0, 1, 2, 3};
    TypeParam test_range{source};

    // pipe notation
    auto v = test_range | seqan3::views::enforce_random_access;
    EXPECT_RANGE_EQ(v, source);

    // function notation
    auto v2 = seqan3::views::enforce_random_access(test_range);
    EXPECT_RANGE_EQ(v2, source);

    // combinability
    auto v3 = test_range | seqan3::views::enforce_random_access | std::views::drop(1);
    EXPECT_RANGE_EQ(v3, (std::vector{1, 2, 3}));
}

// ----------------------------------------------------------------------------
// iterator test
// ----------------------------------------------------------------------------

template <typename rng_t>
struct iterator_fixture<std::type_identity<rng_t>> : public ::testing::Test
{
    static_assert(!std::ranges::random_access_range<rng_t>);

    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = true;

    std::vector<int> expected_range{0, 1, 2, 3, 4, 5, 6, 7};

    rng_t urng{expected_range};
    decltype(urng | seqan3::views::enforce_random_access) test_range = urng | seqan3::views::enforce_random_access;
};

using test_type = ::testing::Types<std::type_identity<common_pseudo_random_access_range>,
                                   std::type_identity<sentinel_pseudo_random_access_range>>;

INSTANTIATE_TYPED_TEST_SUITE_P(pseudo_random_access_view_iterator, iterator_fixture, test_type, );
