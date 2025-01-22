// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alignment/pairwise/edit_distance_unbanded.hpp>

using seqan3::detail::proxy_reference;

template <typename T>
struct remove_reference
{
    using type = std::remove_reference_t<T>;
};

template <typename T>
struct remove_reference<proxy_reference<T>>
{
    using type = T;
};

template <typename T>
using remove_reference_t = typename remove_reference<T>::type;

template <typename reference_t>
constexpr bool is_const_ref_v = std::is_const_v<remove_reference_t<reference_t>>;

template <typename T>
using reference_test = ::testing::Test;

template <typename T>
using proxy_reference_test = ::testing::Test;

using reference_types = ::testing::Types<int &,
                                         int const &,
                                         proxy_reference<int>,
                                         proxy_reference<int> const,
                                         proxy_reference<int const>,
                                         proxy_reference<int const> const>;
using proxy_reference_types = ::testing::Types<proxy_reference<int>, proxy_reference<int const>>;

TYPED_TEST_SUITE(reference_test, reference_types, );
TYPED_TEST_SUITE(proxy_reference_test, proxy_reference_types, );

TYPED_TEST(reference_test, construct_with_lvalue)
{
    using reference_t = TypeParam;
    int a = 5;

    reference_t x{a}; // tracks a

    EXPECT_EQ(a, 5);
    EXPECT_EQ(x, 5); // tracks a

    a = 15;
    EXPECT_EQ(a, 15);
    EXPECT_EQ(x, 15); // tracks a

    if constexpr (!is_const_ref_v<reference_t>)
    {
        x = 115;
        EXPECT_EQ(a, 115);
        EXPECT_EQ(x, 115); // tracks a
    }
}

TYPED_TEST(reference_test, construct_with_reference)
{
    using reference_t = TypeParam;
    int a = 5;

    reference_t x{a}; // tracks a
    reference_t y{x}; // tracks a

    a = 15;
    EXPECT_EQ(a, 15);
    EXPECT_EQ(x, 15); // tracks a
    EXPECT_EQ(y, 15); // tracks a

    if constexpr (!is_const_ref_v<reference_t>)
    {
        x = 115;
        EXPECT_EQ(a, 115);
        EXPECT_EQ(x, 115); // tracks a
        EXPECT_EQ(y, 115); // tracks a

        y = 1115;
        EXPECT_EQ(a, 1115);
        EXPECT_EQ(x, 1115); // tracks a
        EXPECT_EQ(y, 1115); // tracks a
    }
}

TYPED_TEST(reference_test, assign_with_lvalue)
{
    using reference_t = TypeParam;
    int a = 5;

    reference_t x = a; // tracks a

    a = 15;
    EXPECT_EQ(a, 15);
    EXPECT_EQ(x, 15); // tracks a

    if constexpr (!is_const_ref_v<reference_t>)
    {
        x = 115;
        EXPECT_EQ(a, 115);
        EXPECT_EQ(x, 115); // tracks a
    }
}

TYPED_TEST(reference_test, assign_with_reference)
{
    using reference_t = TypeParam;
    int a = 5;

    reference_t x = a; // tracks a
    reference_t y = x; // tracks a

    a = 15;
    EXPECT_EQ(a, 15);
    EXPECT_EQ(x, 15); // tracks a
    EXPECT_EQ(y, 15); // tracks a

    if constexpr (!is_const_ref_v<reference_t>)
    {
        x = 115;
        EXPECT_EQ(a, 115);
        EXPECT_EQ(x, 115); // tracks a
        EXPECT_EQ(y, 115); // tracks a

        y = 1115;
        EXPECT_EQ(a, 1115);
        EXPECT_EQ(x, 1115); // tracks a
        EXPECT_EQ(y, 1115); // tracks a
    }
}

TYPED_TEST(reference_test, assign_value)
{
    using reference_t = TypeParam;
    int a = 5;
    int b = 500;

    reference_t x = a; // tracks a
    reference_t y = b; // tracks b

    a = 15;
    EXPECT_EQ(a, 15);
    EXPECT_EQ(x, 15); // tracks a
    EXPECT_EQ(b, 500);
    EXPECT_EQ(y, 500); // tracks b

    if constexpr (!is_const_ref_v<reference_t>)
    {
        x = y; // still tracks a
        EXPECT_EQ(a, 500);
        EXPECT_EQ(x, 500); // tracks a
        EXPECT_EQ(b, 500);
        EXPECT_EQ(y, 500); // tracks b

        x = 5;
        y = 50;
        EXPECT_EQ(a, 5);
        EXPECT_EQ(x, 5); // tracks a
        EXPECT_EQ(b, 50);
        EXPECT_EQ(y, 50); // tracks b
    }
}

TYPED_TEST(proxy_reference_test, default_construct_and_move_assign)
{
    using reference_t = TypeParam;
    int a = 5;
    int b = 500;

    reference_t x{};    // tracks nothing, yet
    x = reference_t{a}; // tracks a, this statement is not valid for proxy_reference<int> const
    reference_t y = b;  // tracks b

    a = 15;
    EXPECT_EQ(a, 15);
    EXPECT_EQ(x, 15); // tracks a
    EXPECT_EQ(b, 500);
    EXPECT_EQ(y, 500); // tracks b

    x = reference_t{b}; // tracks b, this statement is not valid for proxy_reference<int> const
    EXPECT_EQ(a, 15);
    EXPECT_EQ(x, 500); // tracks b
    EXPECT_EQ(b, 500);
    EXPECT_EQ(y, 500); // tracks b

    if constexpr (!is_const_ref_v<reference_t>)
    {
        x = 5; // has no effect
        y = 50;
        EXPECT_EQ(a, 15);
        EXPECT_EQ(x, 50); // tracks b
        EXPECT_EQ(b, 50);
        EXPECT_EQ(y, 50); // tracks b
    }

    a = 10;
    b = 100;
    EXPECT_EQ(a, 10);
    EXPECT_EQ(x, 100); // tracks b
    EXPECT_EQ(b, 100);
    EXPECT_EQ(y, 100); // tracks b
}

TYPED_TEST(proxy_reference_test, move_construct)
{
    using reference_t = TypeParam;
    int a = 5;
    int b = 500;

    reference_t x{};    // tracks nothing, yet
    x = reference_t{a}; // tracks a, this statement is not valid for proxy_reference<int> const
    reference_t y = b;  // tracks b

    a = 15;
    EXPECT_EQ(a, 15);
    EXPECT_EQ(x, 15); // tracks a
    EXPECT_EQ(b, 500);
    EXPECT_EQ(y, 500); // tracks b

    x = std::move(y); // tracks b, this statement is not valid for proxy_reference<int> const
    EXPECT_EQ(a, 15);
    EXPECT_EQ(x, 500); // tracks b
    EXPECT_EQ(b, 500);
    EXPECT_EQ(y, 500); // tracks b

    if constexpr (!is_const_ref_v<reference_t>)
    {
        x = 5; // has no effect
        y = 50;
        EXPECT_EQ(a, 15);
        EXPECT_EQ(x, 50); // tracks b
        EXPECT_EQ(b, 50);
        EXPECT_EQ(y, 50); // tracks b
    }

    a = 10;
    b = 100;
    EXPECT_EQ(a, 10);
    EXPECT_EQ(x, 100); // tracks b
    EXPECT_EQ(b, 100);
    EXPECT_EQ(y, 100); // tracks b
}
