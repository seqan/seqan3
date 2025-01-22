// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <gtest/gtest.h>

#include <algorithm>
#include <iterator>
#include <ranges>

#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/range/concept.hpp>

struct range_test_fixture
{
    // The value type the std::ranges::range_value_t<range> accepts.
    using range_value_t = void;
    // The value type the std::ranges::range_reference_t<range> accepts.
    using range_reference_t = void;
    // The value type the std::ranges::range_value_t<range const> accepts.
    using range_const_value_t = void;
    // The value type the std::ranges::range_reference_t<range const> accepts.
    using range_const_reference_t = void;

    // Whether the range is a std::ranges::input_range
    static constexpr bool input_range = false;
    // Whether the range is a std::ranges::forward_range
    static constexpr bool forward_range = false;
    // Whether the range is a std::ranges::bidirectional_range
    static constexpr bool bidirectional_range = false;
    // Whether the range is a std::ranges::random_access_range
    static constexpr bool random_access_range = false;
    // Whether the range is a std::ranges::contiguous_range
    static constexpr bool contiguous_range = false;

    // Whether the range is a std::ranges::output_range
    static constexpr bool output_range = false;
    // The value type the std::ranges::output_range accepts.
    using output_value_t = void;

    // Whether the range is a std::ranges::common_range
    static constexpr bool common_range = false;
    // Whether the range is a std::ranges::viewable_range
    static constexpr bool viewable_range = false;
    // Whether the range is a std::ranges::view
    static constexpr bool view = false;
    // Whether the range is a std::ranges::sized_range
    static constexpr bool sized_range = false;
    // Whether the range is a seqan3::const_iterable_range
    static constexpr bool const_iterable_range = false;

    // Whether the range has a size() member.
    static constexpr bool size_member = false;
    // Whether the const range has a size() member.
    static constexpr bool const_size_member = false;
    // Whether the range has a operator[]() member (const version will be checked if range is const_iterable_range).
    static constexpr bool subscript_member = false;

    // The elements the range stores. This should typically be a std::vector.
    auto expected_range() = delete; // to implement
    // The actual range stores.
    auto range() = delete; // to implement

    // This will be called by the subscript_member test, e.g. expect_range_value_equal(range()[0], expected_range()[0])
    template <typename range_value_t, typename expected_range_value_t>
    static void expect_range_value_equal(range_value_t && range_value, expected_range_value_t && expected_range_value)
    {
        static_assert(std::equality_comparable_with<range_value_t, expected_range_value_t>,
                      "The reference types of range_test_fixture::range() and range_test_fixture::expected_range() "
                      "must be equality comparable. If they are not, you may specify a custom void "
                      "expect_range_value_equal(range_value, expected_range_value) function in the fixture.");

        EXPECT_EQ(range_value, expected_range_value);
    }
};

template <typename range_test_fixture_t>
    requires std::derived_from<range_test_fixture_t, range_test_fixture>
struct range_test : public range_test_fixture_t, public ::testing::Test
{};

// ==============================================================

// fwd iterator_fixture it is actually defined in iterator_test_template.hpp
template <typename T>
struct iterator_fixture;

// We can derive iterator_fixture (declared in iterator_test_template.hpp) from this range_test_template and define
// everything s.t. INSTANTIATE_TYPED_TEST_SUITE_P should work without any additional definitions.
template <typename range_test_fixture_t>
    requires std::derived_from<range_test_fixture_t, range_test_fixture>
struct iterator_fixture<range_test_fixture_t> : public ::testing::Test
{
    static constexpr size_t iterator_tag_index = std::max({range_test_fixture_t::input_range ? 0 : 0,
                                                           range_test_fixture_t::forward_range ? 1 : 0,
                                                           range_test_fixture_t::bidirectional_range ? 2 : 0,
                                                           range_test_fixture_t::random_access_range ? 3 : 0,
                                                           range_test_fixture_t::contiguous_range ? 4 : 0});

    using iterator_tag = std::tuple_element_t<iterator_tag_index,
                                              std::tuple<std::input_iterator_tag,
                                                         std::forward_iterator_tag,
                                                         std::bidirectional_iterator_tag,
                                                         std::random_access_iterator_tag,
                                                         std::contiguous_iterator_tag>>;

    static constexpr bool const_iterable = range_test_fixture_t::const_iterable_range;

    template <typename iter_value_t, typename expected_iter_value_t>
    static void expect_eq(iter_value_t && iter_value, expected_iter_value_t && expected_iter_value)
    {
        range_test_fixture_t::expect_range_value_equal(std::forward<iter_value_t>(iter_value),
                                                       std::forward<expected_iter_value_t>(expected_iter_value));
    }

    virtual void SetUp() override
    {
        // re-initialise iterator_fixture after each TestCase in case this is an input iterator
        test_range = range_test_fixture_t{}.range();
    }

    using test_range_t = decltype(range_test_fixture_t{}.range());
    test_range_t test_range;

    using expected_range_t = decltype(range_test_fixture_t{}.expected_range());
    expected_range_t expected_range = range_test_fixture_t{}.expected_range();
};

// ==============================================================

TYPED_TEST_SUITE_P(range_test);

template <typename range_t>
concept has_size_member = requires (range_t range) {
    { range.size() };
};

template <typename range_t>
concept has_subscript_member = requires (range_t range) {
    { range[0] };
};

TYPED_TEST_P(range_test, concept_check)
{
    auto range = this->range();
    using range_t = decltype(range);

    // general range properties
    EXPECT_TRUE(std::ranges::range<range_t>);
    EXPECT_EQ(TestFixture::const_iterable_range, std::ranges::range<range_t const>);
    EXPECT_EQ(TestFixture::const_iterable_range, seqan3::const_iterable_range<range_t>);

    // ranges that have a std::iterator_traits<It>::iterator_concept
    EXPECT_EQ(TestFixture::output_range, (std::ranges::output_range<range_t, typename TestFixture::output_value_t>));

    EXPECT_EQ(TestFixture::input_range, std::ranges::input_range<range_t>);
    EXPECT_EQ(TestFixture::input_range && TestFixture::const_iterable_range, std::ranges::input_range<range_t const>);

    EXPECT_EQ(TestFixture::forward_range, std::ranges::forward_range<range_t>);
    EXPECT_EQ(TestFixture::forward_range && TestFixture::const_iterable_range,
              std::ranges::forward_range<range_t const>);

    EXPECT_EQ(TestFixture::bidirectional_range, std::ranges::bidirectional_range<range_t>);
    EXPECT_EQ(TestFixture::bidirectional_range && TestFixture::const_iterable_range,
              std::ranges::bidirectional_range<range_t const>);

    EXPECT_EQ(TestFixture::random_access_range, std::ranges::random_access_range<range_t>);
    EXPECT_EQ(TestFixture::random_access_range && TestFixture::const_iterable_range,
              std::ranges::random_access_range<range_t const>);

    EXPECT_EQ(TestFixture::contiguous_range, std::ranges::contiguous_range<range_t>);
    EXPECT_EQ(TestFixture::contiguous_range && TestFixture::const_iterable_range,
              std::ranges::contiguous_range<range_t const>);

    // specic range properties that are all orthogonal
    EXPECT_EQ(TestFixture::common_range, std::ranges::common_range<range_t>);
    EXPECT_EQ(TestFixture::common_range && TestFixture::const_iterable_range, std::ranges::common_range<range_t const>);

    EXPECT_EQ(TestFixture::viewable_range, std::ranges::viewable_range<range_t>);
    EXPECT_EQ(TestFixture::viewable_range && TestFixture::const_iterable_range,
              std::ranges::viewable_range<range_t const>);

    EXPECT_EQ(TestFixture::view, std::ranges::view<range_t>);
    // there can't be any view that satisfies std::ranges::view<range_t const>,
    // since it requires std::movable<range_t const> which requires std::assignable_from<range_t const &, range_t const>
    // which is defined as requires(range_t const & lhs, range_t const && rhs)
    // { { lhs = std::forward<range_t const>(rhs) } -> std::same_as<range_t const &>; };
    // (you can't assign anything to a const object)

    EXPECT_EQ(TestFixture::sized_range, std::ranges::sized_range<range_t>);
    EXPECT_EQ(TestFixture::sized_range && TestFixture::const_iterable_range, std::ranges::sized_range<range_t const>);

    // member properties
    EXPECT_EQ(TestFixture::size_member, has_size_member<range_t>);
    EXPECT_EQ(TestFixture::const_size_member, has_size_member<range_t const>);

    EXPECT_EQ(TestFixture::subscript_member, has_subscript_member<range_t>);
    EXPECT_EQ(TestFixture::subscript_member && TestFixture::const_iterable_range, has_subscript_member<range_t const>);

    EXPECT_SAME_TYPE(std::ranges::range_value_t<range_t>, typename TestFixture::range_value_t);
    EXPECT_SAME_TYPE(std::ranges::range_reference_t<range_t>, typename TestFixture::range_reference_t);

    if constexpr (TestFixture::const_iterable_range)
    {
        EXPECT_SAME_TYPE(std::ranges::range_value_t<range_t const>, typename TestFixture::range_const_value_t);
        EXPECT_SAME_TYPE(std::ranges::range_reference_t<range_t const>, typename TestFixture::range_const_reference_t);
    }
}

TYPED_TEST_P(range_test, sized_range)
{
    if constexpr (TestFixture::sized_range)
    {
        {
            auto range = this->range();
            // lvalue
            EXPECT_EQ(std::ranges::size(this->expected_range()), std::ranges::size(range));
            // rvalue
            EXPECT_EQ(std::ranges::size(this->expected_range()), std::ranges::size(this->range()));
        }

        if constexpr (TestFixture::const_iterable_range)
        {
            auto const range = this->range();
            // const lvalue
            EXPECT_EQ(std::ranges::size(this->expected_range()), std::ranges::size(range));
            // const rvalue
            EXPECT_EQ(std::ranges::size(this->expected_range()),
                      std::ranges::size(static_cast<decltype(range) const &&>(this->range())));
        }
    }
}

TYPED_TEST_P(range_test, size_member)
{
    if constexpr (TestFixture::size_member)
    {
        {
            auto range = this->range();
            EXPECT_EQ(std::ranges::size(this->expected_range()), range.size());
        }

        if constexpr (TestFixture::const_iterable_range)
        {
            auto const range = this->range();
            EXPECT_EQ(std::ranges::size(this->expected_range()), range.size());
        }
    }
}

TYPED_TEST_P(range_test, subscript_member)
{
    if constexpr (TestFixture::subscript_member)
    {
        {
            auto range = this->range();
            auto expected_range = this->expected_range();

            for (size_t i = 0; i < std::ranges::size(expected_range); ++i)
            {
                this->expect_range_value_equal(range[i], expected_range[i]);
            }
        }

        if constexpr (TestFixture::const_iterable_range)
        {
            auto const range = this->range();
            auto expected_range = this->expected_range();

            for (size_t i = 0; i < std::ranges::size(expected_range); ++i)
            {
                this->expect_range_value_equal(range[i], expected_range[i]);
            }
        }
    }
}

REGISTER_TYPED_TEST_SUITE_P(range_test, concept_check, sized_range, size_member, subscript_member);
