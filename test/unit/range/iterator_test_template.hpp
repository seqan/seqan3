// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <gtest/gtest.h>

#include <iterator>
#include <memory>
#include <ranges>

#include <seqan3/core/platform.hpp>

template <typename T>
struct iterator_fixture : public ::testing::Test
{
    /* Please provide the following members:
    --------------------------------------------------------------------------------------------------------------------
    using iterator_tag = ...                                   // Defines the iterator functionality you want to test.
                                                               // One of:
                                                               // std::input_iterator_tag
                                                               // std::forward_iterator_tag
                                                               // std::bidirectional_iterator_tag
                                                               // std::random_access_iterator_tag
                                                               // std::contiguous_iterator_tag

    static constexpr bool const_iterable = true/false;         // Also test const_iterability. (const begin/end required)

    t1 test_range;                                             // The range to test the iterators (begin/end required).
    t2 expected_range;                                         // Used to compare the iterator range with.

    --------------------------------------------------------------------------------------------------------------------
    Note: if the reference type of your iterator is not comparable via operator==() to the reference type of
          `expected_range you can additionally specify a custom expect_eq function:

    template <typename A, typename B>
    static void expect_eq(A && begin_iterator_value, B && expected_range_value)
    {
        EXPECT_EQ(begin_iterator_value, expected_range_value);
    }
    */
};

// Helper concept to check whether the test fixture has a member function expect_eq.
template <typename t>
concept has_expect_equal_member_function = requires (t & a) {
    { t::expect_eq(*std::ranges::begin(a.test_range), *std::ranges::begin(a.expected_range)) } -> std::same_as<void>;
};

// Delegates to the test fixture member function `expect_eq` if available and falls back to EXPECT_EQ otherwise.
template <typename T, typename A, typename B>
void expect_iter_value_equal(A && a, B && b)
{
    if constexpr (has_expect_equal_member_function<iterator_fixture<T>>)
        iterator_fixture<std::remove_reference_t<T>>::expect_eq(a, b);
    else
        EXPECT_EQ(a, b);
}

template <typename T, typename it_t, typename rng_it_t>
void expect_iter_equal(it_t && it, rng_it_t && rng_it)
{
    expect_iter_value_equal<T>(*it, *rng_it);
}

// std c++20 input iterator aren't required to have an operator==(iterator_t, iterator_t), but if they have one we
// test the semantic
template <typename type_param_t>
concept iterator_is_equality_comparable =
    std::derived_from<typename iterator_fixture<type_param_t>::iterator_tag, std::forward_iterator_tag>
    || requires (iterator_fixture<type_param_t> & fixture) {
           typename std::ranges::iterator_t<decltype(fixture.test_range)>;

           requires requires (std::ranges::iterator_t<decltype(fixture.test_range)> & it) {
               // we don't assume anything about the return type, this will be done in the tests
               { it == it };
           };
       };

TYPED_TEST_SUITE_P(iterator_fixture);

TYPED_TEST_P(iterator_fixture, concept_check)
{
    using iterator_type = decltype(std::ranges::begin(this->test_range));
    // Ensure that reference types are comparable if no equal_eq function was defined.
    if constexpr (!has_expect_equal_member_function<iterator_fixture<TypeParam>>)
    {
        static_assert(std::equality_comparable_with<decltype(*std::ranges::begin(this->test_range)),
                                                    decltype(*std::ranges::begin(this->expected_range))>,
                      "The reference types of begin_iterator and expected_range must be equality comparable. "
                      "If they are not, you may specify a custom void expect_eq(i1, r2) function in the fixture.");
    }

    // input iterator must always be satisfied
    static_assert(std::input_iterator<decltype(std::ranges::begin(this->expected_range))>,
                  "expected_range must have a begin member function and "
                  "the returned iterator must model std::input_iterator.");
    EXPECT_TRUE(std::input_iterator<iterator_type>);

    EXPECT_EQ(std::forward_iterator<iterator_type>,
              (std::derived_from<typename TestFixture::iterator_tag, std::forward_iterator_tag>));

    EXPECT_EQ(std::bidirectional_iterator<iterator_type>,
              (std::derived_from<typename TestFixture::iterator_tag, std::bidirectional_iterator_tag>));

    EXPECT_EQ(std::random_access_iterator<iterator_type>,
              (std::derived_from<typename TestFixture::iterator_tag, std::random_access_iterator_tag>));

    EXPECT_EQ(std::contiguous_iterator<iterator_type>,
              (std::derived_from<typename TestFixture::iterator_tag, std::contiguous_iterator_tag>));

    if constexpr (TestFixture::const_iterable)
    {
        using const_iterator_type = decltype(std::ranges::cbegin(this->test_range));
        EXPECT_TRUE(std::input_iterator<const_iterator_type>);

        EXPECT_EQ(std::forward_iterator<const_iterator_type>,
                  (std::derived_from<typename TestFixture::iterator_tag, std::forward_iterator_tag>));

        EXPECT_EQ(std::bidirectional_iterator<const_iterator_type>,
                  (std::derived_from<typename TestFixture::iterator_tag, std::bidirectional_iterator_tag>));

        EXPECT_EQ(std::random_access_iterator<const_iterator_type>,
                  (std::derived_from<typename TestFixture::iterator_tag, std::random_access_iterator_tag>));

        EXPECT_EQ(std::contiguous_iterator<const_iterator_type>,
                  (std::derived_from<typename TestFixture::iterator_tag, std::contiguous_iterator_tag>));
    }

    if (!std::derived_from<typename TestFixture::iterator_tag, std::input_iterator_tag>)
    {
        FAIL() << "The iterator tag member type must be one of std::input_iterator_tag, "
               << "std::forward_iterator_tag, std::bidirectional_iterator_tag, std::random_access_iterator_tag, or "
               << "std::contiguous_iterator_tag.";
    }
}

TYPED_TEST_P(iterator_fixture, const_non_const_compatibility)
{
    if constexpr (TestFixture::const_iterable)
    {
        using const_iterator_type = decltype(std::ranges::cbegin(this->test_range));

        [[maybe_unused]] const_iterator_type it{std::ranges::begin(this->test_range)};

        const_iterator_type it2{};
        it2 = std::ranges::begin(this->test_range);

        if constexpr (iterator_is_equality_comparable<TypeParam>)
        {
            EXPECT_EQ(it, it2);
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Input & Forward Iterator
// ---------------------------------------------------------------------------------------------------------------------

TYPED_TEST_P(iterator_fixture, dereference)
{
    expect_iter_equal<TypeParam>(std::ranges::begin(this->test_range), std::ranges::begin(this->expected_range));

    if constexpr (TestFixture::const_iterable)
        expect_iter_equal<TypeParam>(std::ranges::cbegin(this->test_range), std::ranges::begin(this->expected_range));
}

TYPED_TEST_P(iterator_fixture, compare)
{
    EXPECT_FALSE(std::ranges::begin(this->test_range) == std::ranges::end(this->test_range));
    EXPECT_TRUE(std::ranges::begin(this->test_range) != std::ranges::end(this->test_range));
    EXPECT_FALSE(std::ranges::end(this->test_range) == std::ranges::begin(this->test_range));
    EXPECT_TRUE(std::ranges::end(this->test_range) != std::ranges::begin(this->test_range));

    if constexpr (iterator_is_equality_comparable<TypeParam>)
    {
        EXPECT_TRUE(std::ranges::begin(this->test_range) == std::ranges::begin(this->test_range));
        EXPECT_FALSE(std::ranges::begin(this->test_range) != std::ranges::begin(this->test_range));
    }

    if constexpr (TestFixture::const_iterable)
    {
        if constexpr (iterator_is_equality_comparable<TypeParam>)
        {
            EXPECT_TRUE(std::ranges::cbegin(this->test_range) == std::ranges::cbegin(this->test_range));
            EXPECT_FALSE(std::ranges::cbegin(this->test_range) != std::ranges::cbegin(this->test_range));
        }

        EXPECT_FALSE(std::ranges::cbegin(this->test_range) == std::ranges::cend(this->test_range));
        EXPECT_TRUE(std::ranges::cbegin(this->test_range) != std::ranges::cend(this->test_range));
        EXPECT_FALSE(std::ranges::cend(this->test_range) == std::ranges::cbegin(this->test_range));
        EXPECT_TRUE(std::ranges::cend(this->test_range) != std::ranges::cbegin(this->test_range));

        // (non-const lhs)
        if constexpr (iterator_is_equality_comparable<TypeParam>)
        {
            EXPECT_TRUE(std::ranges::begin(this->test_range) == std::ranges::cbegin(this->test_range));
            EXPECT_FALSE(std::ranges::begin(this->test_range) != std::ranges::cbegin(this->test_range));
        }

        EXPECT_FALSE(std::ranges::begin(this->test_range) == std::ranges::cend(this->test_range));
        EXPECT_TRUE(std::ranges::begin(this->test_range) != std::ranges::cend(this->test_range));
        EXPECT_FALSE(std::ranges::end(this->test_range) == std::ranges::cbegin(this->test_range));
        EXPECT_TRUE(std::ranges::end(this->test_range) != std::ranges::cbegin(this->test_range));

        // (non-const rhs)
        if constexpr (iterator_is_equality_comparable<TypeParam>)
        {
            EXPECT_TRUE(std::ranges::cbegin(this->test_range) == std::ranges::begin(this->test_range));
            EXPECT_FALSE(std::ranges::cbegin(this->test_range) != std::ranges::begin(this->test_range));
        }

        EXPECT_FALSE(std::ranges::cend(this->test_range) == std::ranges::begin(this->test_range));
        EXPECT_TRUE(std::ranges::cend(this->test_range) != std::ranges::begin(this->test_range));
        EXPECT_FALSE(std::ranges::cbegin(this->test_range) == std::ranges::end(this->test_range));
        EXPECT_TRUE(std::ranges::cbegin(this->test_range) != std::ranges::end(this->test_range));
    }
}

template <typename test_type, typename it_begin_t, typename it_sentinel_t, typename rng_t>
inline void move_forward_pre_test(it_begin_t && it_begin, it_sentinel_t && it_end, rng_t && rng)
{
    // pre-increment
    auto rng_it = std::ranges::begin(rng);
    auto rng_it_end = std::ranges::end(rng);
    auto it = std::move(it_begin);

    EXPECT_NE(rng_it, rng_it_end);
    EXPECT_NE(it, it_end);

    for (; true;)
    {
        // if it_begin_t is copy_constructible copy result, otherwise take it by reference (if move-only iterator)
        using it_copy_or_reference_t =
            std::conditional_t<std::copy_constructible<it_begin_t>, std::remove_reference_t<it_begin_t>, it_begin_t &>;
        it_copy_or_reference_t it_copy_or_reference = ++it;
        ++rng_it;

        if (it == it_end || rng_it == rng_it_end)
            break;

        expect_iter_equal<test_type>(it_copy_or_reference, rng_it);
    }
    EXPECT_EQ(rng_it, rng_it_end);
    EXPECT_EQ(it, it_end);
}

template <typename test_type, typename it_begin_t, typename it_sentinel_t, typename rng_t>
inline void move_forward_post_test(it_begin_t && it_begin, it_sentinel_t && it_end, rng_t && rng)
{
    // post-increment
    auto rng_it = std::ranges::begin(rng);
    auto rng_it_end = std::ranges::end(rng);
    auto it = std::move(it_begin);

    EXPECT_NE(rng_it, rng_it_end);
    EXPECT_NE(it, it_end);

    static constexpr bool is_cpp20_input_iterator = std::same_as<decltype(it_begin++), void>;

    if constexpr (is_cpp20_input_iterator)
    {
        // input iterator can return void for post-increment (expressed by std::weakly_incrementable)
        EXPECT_TRUE(std::input_iterator<it_begin_t>);
        // forward iterator require std::incrementable which requires `{ i++ } -> same_as<I>;`
        EXPECT_FALSE(std::forward_iterator<it_begin_t>);
    }

    for (; it != it_end && rng_it != rng_it_end;)
    {
        expect_iter_equal<test_type>(it, rng_it);

        if constexpr (!is_cpp20_input_iterator)
        {
            expect_iter_equal<test_type>(it++, rng_it++);
        }
        else
        {
            it++;
            rng_it++;
        }
    }
    EXPECT_EQ(rng_it, rng_it_end);
    EXPECT_EQ(it, it_end);
}

TYPED_TEST_P(iterator_fixture, move_forward_pre)
{
    move_forward_pre_test<TypeParam>(std::ranges::begin(this->test_range),
                                     std::ranges::end(this->test_range),
                                     this->expected_range);

    // iterate over it again
    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::forward_iterator_tag>)
        move_forward_pre_test<TypeParam>(std::ranges::begin(this->test_range),
                                         std::ranges::end(this->test_range),
                                         this->expected_range);
}

TYPED_TEST_P(iterator_fixture, move_forward_pre_const)
{
    if constexpr (TestFixture::const_iterable)
    {
        move_forward_pre_test<TypeParam>(std::ranges::cbegin(this->test_range),
                                         std::ranges::cend(this->test_range),
                                         this->expected_range);

        // iterate over it again
        if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::forward_iterator_tag>)
            move_forward_pre_test<TypeParam>(std::ranges::cbegin(this->test_range),
                                             std::ranges::cend(this->test_range),
                                             this->expected_range);
    }
}

TYPED_TEST_P(iterator_fixture, move_forward_post)
{
    move_forward_post_test<TypeParam>(std::ranges::begin(this->test_range),
                                      std::ranges::end(this->test_range),
                                      this->expected_range);

    // iterate over it again
    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::forward_iterator_tag>)
        move_forward_post_test<TypeParam>(std::ranges::begin(this->test_range),
                                          std::ranges::end(this->test_range),
                                          this->expected_range);
}

TYPED_TEST_P(iterator_fixture, move_forward_post_const)
{
    if constexpr (TestFixture::const_iterable)
    {
        move_forward_post_test<TypeParam>(std::ranges::cbegin(this->test_range),
                                          std::ranges::cend(this->test_range),
                                          this->expected_range);

        // iterate over it again
        if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::forward_iterator_tag>)
            move_forward_post_test<TypeParam>(std::ranges::cbegin(this->test_range),
                                              std::ranges::cend(this->test_range),
                                              this->expected_range);
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Bidirectional Iterator
// ---------------------------------------------------------------------------------------------------------------------

template <typename it_begin_t, typename it_sentinel_t, typename rng_t>
inline auto last_iterators(it_begin_t const & it_begin, it_sentinel_t const & it_end, rng_t && rng)
{
    it_begin_t it = it_begin;
    auto rng_it = std::ranges::begin(rng);

    for (auto const rng_it_end = std::ranges::end(rng);
         std::ranges::next(it) != it_end && std::ranges::next(rng_it) != rng_it_end;
         ++it, ++rng_it)
        ;

    return std::pair{it, rng_it};
}

template <typename test_type, typename it_begin_t, typename it_sentinel_t, typename rng_t>
inline void move_backward_pre_test(it_begin_t && it_begin, it_sentinel_t && it_end, rng_t && rng)
{
    // move to last position
    auto && [last_it, rng_last_it] = last_iterators(it_begin, it_end, rng);
    auto const rng_it_begin = std::ranges::begin(rng);

    // pre-decrement
    auto it = last_it;
    auto rng_it = rng_last_it;
    for (; it != it_begin && rng_it != rng_it_begin; --rng_it)
    {
        expect_iter_equal<test_type>(it, rng_it);
        --it;
    }

    expect_iter_equal<test_type>(it_begin, rng_it_begin);
}

template <typename test_type, typename it_begin_t, typename it_sentinel_t, typename rng_t>
inline void move_backward_post_test(it_begin_t && it_begin, it_sentinel_t && it_end, rng_t && rng)
{
    // move to last position
    auto && [last_it, rng_last_it] = last_iterators(it_begin, it_end, rng);
    auto const rng_it_begin = std::ranges::begin(rng);

    // post-decrement
    auto it = last_it;
    auto rng_it = rng_last_it;
    for (; it != it_begin && rng_it != rng_it_begin; --rng_it)
    {
        expect_iter_equal<test_type>(it--, rng_it);
    }

    expect_iter_equal<test_type>(it_begin, rng_it_begin);
}

TYPED_TEST_P(iterator_fixture, move_backward_pre)
{
    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::bidirectional_iterator_tag>)
    {
        move_backward_pre_test<TypeParam>(std::ranges::begin(this->test_range),
                                          std::ranges::end(this->test_range),
                                          this->expected_range);

        if constexpr (TestFixture::const_iterable)
            move_backward_pre_test<TypeParam>(std::ranges::cbegin(this->test_range),
                                              std::ranges::cend(this->test_range),
                                              this->expected_range);
    }
}

TYPED_TEST_P(iterator_fixture, move_backward_post)
{
    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::bidirectional_iterator_tag>)
    {
        move_backward_post_test<TypeParam>(std::ranges::begin(this->test_range),
                                           std::ranges::end(this->test_range),
                                           this->expected_range);

        if constexpr (TestFixture::const_iterable)
            move_backward_post_test<TypeParam>(std::ranges::cbegin(this->test_range),
                                               std::ranges::cend(this->test_range),
                                               this->expected_range);
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Random Access Iterator
// ---------------------------------------------------------------------------------------------------------------------

template <typename test_type, typename it_begin_t, typename rng_t>
inline void jump_forward_test(it_begin_t && it_begin, rng_t && rng)
{
    size_t sz = std::ranges::distance(rng);
    auto rng_it_begin = std::ranges::begin(rng);

    // Forward
    for (size_t n = 0; n < sz; ++n)
    {
        auto it = it_begin;
        expect_iter_equal<test_type>(it += n, rng_it_begin + n);
        expect_iter_equal<test_type>(it, rng_it_begin + n);
    }

    // Forward copy
    for (size_t n = 0; n < sz; ++n)
    {
        expect_iter_equal<test_type>(it_begin + n, rng_it_begin + n);
        expect_iter_equal<test_type>(it_begin, rng_it_begin);
    }

    // Forward copy friend
    for (size_t n = 0; n < sz; ++n)
    {
        expect_iter_equal<test_type>(n + it_begin, rng_it_begin + n);
        expect_iter_equal<test_type>(it_begin, rng_it_begin);
    }
}

TYPED_TEST_P(iterator_fixture, jump_forward)
{
    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        jump_forward_test<TypeParam>(std::ranges::begin(this->test_range), this->expected_range);

        if constexpr (TestFixture::const_iterable)
            jump_forward_test<TypeParam>(std::ranges::cbegin(this->test_range), std::as_const(this->expected_range));
    }
}

template <typename test_type, typename it_begin_t, typename rng_t>
inline void jump_backward_test(it_begin_t && it_begin, rng_t && rng)
{
    size_t sz = std::ranges::distance(rng);
    auto rng_it_begin = std::ranges::begin(rng);

    auto pre_end_it = it_begin + (sz - 1);
    auto pre_end_rng_it = rng_it_begin + (sz - 1);

    // Backward
    for (size_t n = 0; n < sz; ++n)
    {
        auto it = pre_end_it;
        expect_iter_equal<test_type>(it -= n, pre_end_rng_it - n);
        expect_iter_equal<test_type>(it, pre_end_rng_it - n);
    }

    // Backward copy
    for (size_t n = 0; n < sz; ++n)
    {
        expect_iter_equal<test_type>(pre_end_it - n, pre_end_rng_it - n);
        expect_iter_equal<test_type>(pre_end_it, pre_end_rng_it);
    }

    // Backward copy it + (-n)
    for (size_t n = 0; n < sz; ++n)
    {
        expect_iter_equal<test_type>(pre_end_it + (-1 * n), pre_end_rng_it - n);
        expect_iter_equal<test_type>(pre_end_it, pre_end_rng_it);
    }

    // Backward copy friend through (-n) + it
    for (size_t n = 0; n < sz; ++n)
    {
        expect_iter_equal<test_type>((-1 * n) + pre_end_it, pre_end_rng_it - n);
        expect_iter_equal<test_type>(pre_end_it, pre_end_rng_it);
    }
}

TYPED_TEST_P(iterator_fixture, jump_backward)
{
    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        jump_backward_test<TypeParam>(std::ranges::begin(this->test_range), this->expected_range);

        if constexpr (TestFixture::const_iterable)
            jump_backward_test<TypeParam>(std::ranges::cbegin(this->test_range), std::as_const(this->expected_range));
    }
}

template <typename test_type, typename it_begin_t, typename rng_t>
inline void jump_random_test(it_begin_t && it_begin, rng_t && rng)
{
    size_t sz = std::ranges::distance(rng);

    for (size_t n = 0; n < sz; ++n)
        expect_iter_value_equal<test_type>(it_begin[n], rng[n]);
}

TYPED_TEST_P(iterator_fixture, jump_random)
{
    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        jump_random_test<TypeParam>(std::ranges::begin(this->test_range), this->expected_range);

        if constexpr (TestFixture::const_iterable)
            jump_random_test<TypeParam>(std::ranges::cbegin(this->test_range), std::as_const(this->expected_range));
    }
}

template <typename iterator_t, typename rng_t>
inline void difference_test(iterator_t && it_begin, iterator_t && it_end, rng_t && rng)
{
    using difference_t = std::iter_difference_t<iterator_t>;
    difference_t size = std::ranges::distance(rng);

    for (difference_t n = 0; n <= size; ++n)
    {
        EXPECT_EQ(n, (it_begin + n) - it_begin);
        EXPECT_EQ(-n, it_begin - (it_begin + n));
    }

    for (difference_t n = 0; n <= size; ++n)
    {
        EXPECT_EQ(n, it_end - (it_end - n));
        EXPECT_EQ(-n, (it_end - n) - it_end);
    }
}

TYPED_TEST_P(iterator_fixture, difference_common)
{
    static constexpr bool is_random_access =
        std::derived_from<typename TestFixture::iterator_tag, std::random_access_iterator_tag>;

    if constexpr (is_random_access)
    {
        auto it = std::ranges::begin(this->test_range);
        auto sentinel = std::ranges::next(it, std::ranges::end(this->test_range));
        difference_test(it, sentinel, this->expected_range);
    }

    if constexpr (is_random_access && TestFixture::const_iterable)
    {
        auto const_it = std::ranges::cbegin(this->test_range);
        auto const_sentinel = std::ranges::next(const_it, std::ranges::cend(this->test_range));
        difference_test(const_it, const_sentinel, std::as_const(this->expected_range));
    }
}

TYPED_TEST_P(iterator_fixture, difference_sentinel)
{
    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        using difference_t = std::ranges::range_difference_t<decltype(this->test_range)>;

        auto && begin = std::ranges::begin(this->test_range);
        auto && end = std::ranges::end(this->test_range);
        difference_t size = std::ranges::distance(this->expected_range);

        EXPECT_EQ(size, end - begin);
        EXPECT_EQ(-size, begin - end);

        if constexpr (TestFixture::const_iterable)
        {
            auto && cbegin = std::ranges::cbegin(this->test_range);
            auto && cend = std::ranges::cend(this->test_range);

            EXPECT_EQ(size, cend - cbegin);
            EXPECT_EQ(-size, cbegin - cend);

            EXPECT_EQ(size, end - cbegin);
            EXPECT_EQ(-size, cbegin - end);

            EXPECT_EQ(size, cend - begin);
            EXPECT_EQ(-size, begin - cend);
        }
    }
}

TYPED_TEST_P(iterator_fixture, compare_less)
{
    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        EXPECT_FALSE(std::ranges::begin(this->test_range) < std::ranges::begin(this->test_range));
        EXPECT_TRUE(std::ranges::begin(this->test_range) < std::ranges::next(std::ranges::begin(this->test_range)));
    }

    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::random_access_iterator_tag>
                  && TestFixture::const_iterable)
    {
        EXPECT_FALSE(std::ranges::cbegin(this->test_range) < std::ranges::cbegin(this->test_range));
        EXPECT_TRUE(std::ranges::cbegin(this->test_range) < std::ranges::next(std::ranges::cbegin(this->test_range)));

        // mix
        EXPECT_FALSE(std::ranges::begin(this->test_range) < std::ranges::cbegin(this->test_range));
        EXPECT_TRUE(std::ranges::begin(this->test_range) < std::ranges::next(std::ranges::cbegin(this->test_range)));
        EXPECT_FALSE(std::ranges::cbegin(this->test_range) < std::ranges::begin(this->test_range));
        EXPECT_TRUE(std::ranges::cbegin(this->test_range) < std::ranges::next(std::ranges::begin(this->test_range)));
    }
}

TYPED_TEST_P(iterator_fixture, compare_greater)
{
    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        EXPECT_FALSE(std::ranges::begin(this->test_range) > std::ranges::begin(this->test_range));
        EXPECT_FALSE(std::ranges::begin(this->test_range) > std::ranges::next(std::ranges::begin(this->test_range)));
    }

    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::random_access_iterator_tag>
                  && TestFixture::const_iterable)
    {
        EXPECT_FALSE(std::ranges::cbegin(this->test_range) > std::ranges::cbegin(this->test_range));
        EXPECT_FALSE(std::ranges::cbegin(this->test_range) > std::ranges::next(std::ranges::cbegin(this->test_range)));

        // mix
        EXPECT_FALSE(std::ranges::begin(this->test_range) > std::ranges::cbegin(this->test_range));
        EXPECT_FALSE(std::ranges::begin(this->test_range) > std::ranges::next(std::ranges::cbegin(this->test_range)));
        EXPECT_FALSE(std::ranges::cbegin(this->test_range) > std::ranges::begin(this->test_range));
        EXPECT_FALSE(std::ranges::cbegin(this->test_range) > std::ranges::next(std::ranges::begin(this->test_range)));
    }
}

TYPED_TEST_P(iterator_fixture, compare_leq)
{
    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        EXPECT_TRUE(std::ranges::begin(this->test_range) <= std::ranges::begin(this->test_range));
        EXPECT_TRUE(std::ranges::begin(this->test_range) <= std::ranges::next(std::ranges::begin(this->test_range)));
    }

    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::random_access_iterator_tag>
                  && TestFixture::const_iterable)
    {
        EXPECT_TRUE(std::ranges::cbegin(this->test_range) <= std::ranges::cbegin(this->test_range));
        EXPECT_TRUE(std::ranges::cbegin(this->test_range) <= std::ranges::next(std::ranges::cbegin(this->test_range)));

        // mix
        EXPECT_TRUE(std::ranges::begin(this->test_range) <= std::ranges::cbegin(this->test_range));
        EXPECT_TRUE(std::ranges::begin(this->test_range) <= std::ranges::next(std::ranges::cbegin(this->test_range)));
        EXPECT_TRUE(std::ranges::cbegin(this->test_range) <= std::ranges::begin(this->test_range));
        EXPECT_TRUE(std::ranges::cbegin(this->test_range) <= std::ranges::next(std::ranges::begin(this->test_range)));
    }
}

TYPED_TEST_P(iterator_fixture, compare_geq)
{
    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        EXPECT_TRUE(std::ranges::begin(this->test_range) >= std::ranges::begin(this->test_range));
        EXPECT_FALSE(std::ranges::begin(this->test_range) >= std::ranges::next(std::ranges::begin(this->test_range)));
    }

    if constexpr (std::derived_from<typename TestFixture::iterator_tag, std::random_access_iterator_tag>
                  && TestFixture::const_iterable)
    {
        EXPECT_TRUE(std::ranges::cbegin(this->test_range) >= std::ranges::cbegin(this->test_range));
        EXPECT_FALSE(std::ranges::cbegin(this->test_range) >= std::ranges::next(std::ranges::cbegin(this->test_range)));

        // mix
        EXPECT_TRUE(std::ranges::begin(this->test_range) >= std::ranges::cbegin(this->test_range));
        EXPECT_FALSE(std::ranges::begin(this->test_range) >= std::ranges::next(std::ranges::cbegin(this->test_range)));
        EXPECT_TRUE(std::ranges::cbegin(this->test_range) >= std::ranges::begin(this->test_range));
        EXPECT_FALSE(std::ranges::cbegin(this->test_range) >= std::ranges::next(std::ranges::begin(this->test_range)));
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Contiguous Iterator
// ---------------------------------------------------------------------------------------------------------------------

template <typename iterator_t>
inline void address_difference_test(iterator_t it_begin, iterator_t it_end)
{
    // contiguous_iterator only requires to_address of the iterator_t, but not sentinel_t.
    using difference_t = std::iter_difference_t<iterator_t>;

    difference_t const size = it_end - it_begin;

    for (difference_t i = 0u; i <= size; ++i)
    {
        iterator_t it = it_begin + i;

        if (it != it_end)
        {
            // https://eel.is/c++draft/iterator.concept.contiguous#2.1
            // to_address(a) == addressof(*a)
            EXPECT_EQ(std::to_address(it), std::addressof(*it));
        }

        // https://eel.is/c++draft/iterator.concept.contiguous#2.2
        // to_address(b) == to_address(a) + D(b - a)
        // to_address(c) == to_address(a) + D(c - a)
        EXPECT_EQ(std::to_address(it), std::to_address(it_begin) + static_cast<difference_t>(it - it_begin));
        EXPECT_EQ(std::to_address(it), std::to_address(it_end) - static_cast<difference_t>(it_end - it));
    }
}

TYPED_TEST_P(iterator_fixture, address_difference)
{
    static constexpr bool is_contiguous =
        std::derived_from<typename TestFixture::iterator_tag, std::contiguous_iterator_tag>;

    if constexpr (is_contiguous)
    {
        auto it = std::ranges::begin(this->test_range);
        auto sentinel_it = std::ranges::next(it, std::ranges::end(this->test_range));
        address_difference_test(it, sentinel_it);
    }

    if constexpr (is_contiguous && TestFixture::const_iterable)
    {
        auto it = std::ranges::cbegin(this->test_range);
        auto sentinel_it = std::ranges::next(it, std::ranges::cend(this->test_range));
        address_difference_test(it, sentinel_it);
    }
}

REGISTER_TYPED_TEST_SUITE_P(iterator_fixture,
                            concept_check,
                            const_non_const_compatibility,
                            dereference,
                            compare,
                            move_forward_pre,
                            move_forward_pre_const,
                            move_forward_post,
                            move_forward_post_const,
                            move_backward_pre,
                            move_backward_post,
                            jump_forward,
                            jump_backward,
                            jump_random,
                            difference_common,
                            difference_sentinel,
                            compare_less,
                            compare_greater,
                            compare_leq,
                            compare_geq,
                            address_difference);
