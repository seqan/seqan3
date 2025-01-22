// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <deque>
#include <forward_list>
#include <list>
#include <ranges>
#include <vector>

#include <seqan3/utility/range/to.hpp>

/* These test cases are taken from the proposal https://wg21.link/p1206r7.
 * This demonstrates the subset of syntaxes constructs that compile.
 */
TEST(to_test, overview)
{
    auto l = std::views::iota(1, 10);

    // create a vector with the elements of l
    auto vec = seqan3::ranges::to<std::vector<int>>(l);

    // Specify an allocator
    auto b = seqan3::ranges::to<std::vector<int, std::allocator<int>>>(l, std::allocator<int>{});

    // deducing value_type
    auto c = seqan3::ranges::to<std::vector>(l);

    // explicit conversion int -> long
    auto d = seqan3::ranges::to<std::vector<long>>(l);

    // Supports converting associative container to sequence containers
    //auto f = seqan3::ranges::to<std::vector>(m); // Not supported by seqan3

    // Supports converting sequence containers to associative ones
    //auto g = seqan3::ranges::to<std::map>(f); // Not supported by seqan3

    // Pipe syntaxe
    auto g = l | std::ranges::views::take(42) | seqan3::ranges::to<std::vector>();

    // Pipe syntax with allocator
    auto h = l | std::ranges::views::take(42) | seqan3::ranges::to<std::vector>(std::allocator<int>{});

    // The pipe syntax also support specifying the type and conversions
    auto i = l | std::ranges::views::take(42) | seqan3::ranges::to<std::vector<long>>();

    // Nested ranges
    std::list<std::forward_list<int>> lst = {{0, 1, 2, 3}, {4, 5, 6, 7}};
    auto vec1 = seqan3::ranges::to<std::vector<std::vector<int>>>(lst);
    auto vec2 = seqan3::ranges::to<std::vector<std::deque<double>>>(lst);
}

// Check that converting a range to `std::vector<int>` works with the function call syntax.
TEST(to_test, function_call_explicit_vector)
{
    auto lst = std::views::iota(1, 10);
    auto vec = seqan3::ranges::to<std::vector<int>>(lst);
    auto expected = std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9};
    EXPECT_EQ(vec, expected);
}

/* Check that converting a range to `std::vector<int, UserDefinedAllocator>` works with the function call
 * syntax and an user-defined allocator.
 */
TEST(to_test, function_call_explicit_vector_with_allocator)
{
    auto lst = std::views::iota(1, 10);
    auto vec = seqan3::ranges::to<std::vector<int, std::allocator<int>>>(lst, std::allocator<int>{});
    auto expected = std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9};
    EXPECT_EQ(vec, expected);
}

// Check that converting a range to `std::vector` works with the function call syntax.
TEST(to_test, function_call_implicit_vector)
{
    auto lst = std::views::iota(1, 10);
    auto vec = seqan3::ranges::to<std::vector>(lst);
    auto expected = std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9};
    EXPECT_EQ(vec, expected);
}

// Check that converting a range to `std::vector` works with the function call syntax and an user-defined allocator.
TEST(to_test, function_call_implicit_vector_with_allocator)
{
    auto lst = std::views::iota(1, 10);
    auto vec = seqan3::ranges::to<std::vector>(lst, std::allocator<int>{});
    auto expected = std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9};
    EXPECT_EQ(vec, expected);
}

// Check that converting a range to `std::vector<int>` works with the function call syntax.
TEST(to_test, function_call_explicit_vector_with_conversion)
{
    auto lst = std::views::iota(1, 10);
    auto vec = seqan3::ranges::to<std::vector<double>>(lst);
    auto expected = std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9};
    EXPECT_EQ(vec, expected);
}

// Check that converting a range to `std::vector<int>` works using pipe syntax.
TEST(to_test, pipe_syntax_explicit_vector)
{
    auto lst = std::views::iota(1, 10);
    auto vec = lst | seqan3::ranges::to<std::vector<int>>();
    auto expected = std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9};
    EXPECT_EQ(vec, expected);
}

/* Check that converting a range to `std::vector<int, UserDefinedAllocator>` works with pipe
 * syntax and an user-defined allocator.
 */
TEST(to_test, pipe_syntax_explicit_vector_with_allocator)
{
    auto lst = std::views::iota(1, 10);
    auto vec = lst | seqan3::ranges::to<std::vector<int, std::allocator<int>>>(std::allocator<int>{});
    auto expected = std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9};
    EXPECT_EQ(vec, expected);
}

// Check that converting a range to `std::vector` works using pipe syntax.
TEST(to_test, pipe_syntax_implicit_vector)
{
    auto lst = std::views::iota(1, 10);
    auto vec = lst | seqan3::ranges::to<std::vector>();
    auto expected = std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9};
    EXPECT_EQ(vec, expected);
}

// Check that converting a range to `std::vector` works with pipe syntax and an user-defined allocator.
TEST(to_test, pipe_syntax_implicit_vector_with_allocator)
{
    auto lst = std::views::iota(1, 10);
    auto vec = lst | seqan3::ranges::to<std::vector>(std::allocator<int>{});
    auto expected = std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9};
    EXPECT_EQ(vec, expected);
}

// Check that converting a range to `std::vector<int>` works using pipe syntax.
TEST(to_test, pipe_syntax_explicit_vector_with_conversion)
{
    auto lst = std::views::iota(1, 10);
    auto vec = lst | seqan3::ranges::to<std::vector<double>>();
    auto expected = std::vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9};
    EXPECT_EQ(vec, expected);
}

// Check that converting nested ranges can be converted to nested containers using function call syntax.
TEST(to_test, syntax_explicit_vector)
{
    auto lst = std::list<std::forward_list<int>>{{1, 2, 3}, {4, 5, 6, 7}};
    auto vec = seqan3::ranges::to<std::vector<std::vector<int>>>(lst);
    auto expected = std::vector<std::vector<int>>{{1, 2, 3}, {4, 5, 6, 7}};
    EXPECT_EQ(vec, expected);
}

// Check other converting types.
TEST(to_test, various_types)
{
    auto lst = std::views::iota(1, 10);
    EXPECT_EQ(lst | seqan3::ranges::to<std::vector<int>>(), (std::vector<int>{1, 2, 3, 4, 5, 6, 7, 8, 9}));
    EXPECT_EQ(lst | seqan3::ranges::to<std::list<int>>(), (std::list<int>{1, 2, 3, 4, 5, 6, 7, 8, 9}));
    EXPECT_EQ(lst | seqan3::ranges::to<std::deque<int>>(), (std::deque<int>{1, 2, 3, 4, 5, 6, 7, 8, 9}));
}
