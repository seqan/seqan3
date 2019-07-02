// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <gtest/gtest.h>

#include <seqan3/std/iterator>

using namespace seqan3;

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

    static constexpr bool const_iterable = true/false;         // Also test const_iterability.

    range_t range_to_compare;                                  // Used to compare the iterator range with.

    T begin_iterator;                                          // Begin of the iterator (range) to test.
    T end_iterator;                                            // End of the iterator (range) to test.

    --------------------------------------------------------------------------------------------------------------------
    // Only need to be specified if const_iterable == true.
    --------------------------------------------------------------------------------------------------------------------
    const_T begin_const_iterator;                               // Begin of the const iterator (range) to test.
    const_T end_const_iterator;                                 // end of the const iterator (range) to test.

    --------------------------------------------------------------------------------------------------------------------
    Note: if the reference type of your iterator is not comparable via operator==() to the reference type of
          `range_to_compare you can additionally specify a custom expect_eq function:

    template <typename A, typename B>
    static void expect_eq(A && begin_iterator_value, B && range_to_compare_value)
    {
        EXPECT_EQ(begin_iterator_value, range_to_compare_value);
    }
    */
};

// Helper concept to check whether the test fixture has a member function expect_eq.
template <typename T>
SEQAN3_CONCEPT HasExpectEqualMemberFunction = requires(T a) {
    { a.expect_eq(*a.begin_iterator, *a.range_to_compare.begin()) } -> void;
};

// Delegates to the test fixture member function `expect_eq` if available and falls back to EXPECT_EQ otherwise.
template <typename T, typename A, typename B>
void expext_eq(A && a, B && b)
{
    if constexpr (HasExpectEqualMemberFunction<iterator_fixture<T>>)
        iterator_fixture<std::remove_reference_t<T>>::expect_eq(a, b);
    else
        EXPECT_EQ(a, b);
}

TYPED_TEST_CASE_P(iterator_fixture);

TYPED_TEST_P(iterator_fixture, concept_check)
{
    // Ensure that reference types are comparable if no equal_eq function was defined.
    if constexpr (!HasExpectEqualMemberFunction<iterator_fixture<TypeParam>>)
    {
        static_assert(std::EqualityComparableWith<decltype(*this->begin_iterator),
                                                  decltype(*this->range_to_compare.begin())>,
                      "The reference types of begin_iterator and range_to_compare must be equality comparable. "
                      "If they are not, you may specify a custom void expect_eq(i1, r2) function in the fixture.");
    }

    if constexpr (std::Same<typename TestFixture::iterator_tag, std::input_iterator_tag>)
    {
        static_assert(std::InputIterator<decltype(this->range_to_compare.begin())>,
                      "range_to_compare must have a begin member function and "
                      "the returned iterator must model std::InputIterator.");
        EXPECT_TRUE(std::InputIterator<TypeParam>);
    }
    else if constexpr (std::Same<typename TestFixture::iterator_tag, std::forward_iterator_tag>)
    {
        static_assert(std::ForwardIterator<decltype(this->range_to_compare.begin())>,
                      "range_to_compare must have a begin member function and "
                      "the returned iterator must model std::ForwardIterator.");
        EXPECT_TRUE(std::ForwardIterator<TypeParam>);

        if constexpr (TestFixture::const_iterable)
        {
            EXPECT_TRUE(std::ForwardIterator<decltype(this->begin_const_iterator)>);
        }
    }
    else if constexpr (std::Same<typename TestFixture::iterator_tag, std::bidirectional_iterator_tag>)
    {
        static_assert(std::BidirectionalIterator<decltype(this->range_to_compare.begin())>,
                      "range_to_compare must have a begin member function and "
                      "the returned iterator must model std::BidirectionalIterator.");
        EXPECT_TRUE(std::BidirectionalIterator<TypeParam>);

        if constexpr (TestFixture::const_iterable)
        {
            EXPECT_TRUE(std::BidirectionalIterator<decltype(this->begin_const_iterator)>);
        }
    }
    else if constexpr (std::Same<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        static_assert(std::RandomAccessIterator<decltype(this->range_to_compare.begin())>,
                      "range_to_compare must have a begin member function and "
                      "the returned iterator must model std::RandomAccessIterator.");
        EXPECT_TRUE(std::RandomAccessIterator<TypeParam>);

        if constexpr (TestFixture::const_iterable)
        {
            EXPECT_TRUE(std::RandomAccessIterator<decltype(this->begin_const_iterator)>);
        }
    }
    else
    {
        FAIL() << "The iterator tag member type must be one of std::forward_iterator_tag,"
               << "std::bidirectional_iterator_tag, std::random_access_iterator_tag.";
    }
}

TYPED_TEST_P(iterator_fixture, const_non_const_compatibility)
{
    if constexpr (TestFixture::const_iterable)
    {
        using const_iterator_type = decltype(this->begin_const_iterator);

        const_iterator_type it{this->begin_iterator};

        const_iterator_type it2{};
        it2 = this->begin_iterator;

        EXPECT_EQ(it, it2);
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Input & Forward Iterator
// ---------------------------------------------------------------------------------------------------------------------

TYPED_TEST_P(iterator_fixture, dereference)
{
    expext_eq<TypeParam>(*this->begin_iterator, *this->range_to_compare.begin());

    if constexpr (TestFixture::const_iterable)
        expext_eq<TypeParam>(*this->begin_const_iterator, *this->range_to_compare.begin());
}

TYPED_TEST_P(iterator_fixture, compare)
{
    EXPECT_FALSE(this->begin_iterator == this->end_iterator);
    EXPECT_TRUE(this->begin_iterator  != this->end_iterator);
    EXPECT_FALSE(this->end_iterator   == this->begin_iterator);
    EXPECT_TRUE(this->end_iterator    != this->begin_iterator);

    if constexpr (std::Same<typename TestFixture::iterator_tag, std::forward_iterator_tag>) // iterate over it again
    {
        EXPECT_TRUE(this->begin_iterator  == this->begin_iterator);
        EXPECT_FALSE(this->begin_iterator != this->begin_iterator);
    }

    if constexpr (TestFixture::const_iterable)
    {
        EXPECT_TRUE(this->begin_const_iterator  == this->begin_const_iterator);
        EXPECT_FALSE(this->begin_const_iterator != this->begin_const_iterator);
        EXPECT_FALSE(this->begin_const_iterator == this->end_const_iterator);
        EXPECT_TRUE(this->begin_const_iterator  != this->end_const_iterator);
        EXPECT_FALSE(this->end_const_iterator   == this->begin_const_iterator);
        EXPECT_TRUE(this->end_const_iterator    != this->begin_const_iterator);

        // (non-const lhs)
        EXPECT_TRUE(this->begin_iterator  == this->begin_const_iterator);
        EXPECT_FALSE(this->begin_iterator != this->begin_const_iterator);
        EXPECT_FALSE(this->begin_iterator == this->end_const_iterator);
        EXPECT_TRUE(this->begin_iterator  != this->end_const_iterator);
        EXPECT_FALSE(this->end_iterator   == this->begin_const_iterator);
        EXPECT_TRUE(this->end_iterator    != this->begin_const_iterator);

        // (non-const rhs)
        EXPECT_TRUE(this->begin_const_iterator  == this->begin_iterator);
        EXPECT_FALSE(this->begin_const_iterator != this->begin_iterator);
        EXPECT_FALSE(this->end_const_iterator   == this->begin_iterator);
        EXPECT_TRUE(this->end_const_iterator    != this->begin_iterator);
        EXPECT_FALSE(this->begin_const_iterator == this->end_iterator);
        EXPECT_TRUE(this->begin_const_iterator  != this->end_iterator);
    }
}

template <typename test_type, typename it_begin_t, typename it_sentinel_t, typename rng_t>
inline void move_forward_pre_test(it_begin_t && it_begin, it_sentinel_t && it_end, rng_t && rng)
{
    // pre-increment
    auto rng_it = rng.begin();
    for (auto it = it_begin; it != it_end; ++it, ++rng_it)
        expext_eq<test_type>(*it, *rng_it);
}

template <typename test_type, typename it_begin_t, typename it_sentinel_t, typename rng_t>
inline void move_forward_post_test(it_begin_t && it_begin, it_sentinel_t && it_end, rng_t && rng)
{
    // post-increment
    auto rng_it = rng.begin();
    for (auto it = it_begin; it != it_end; it++, ++rng_it)
        expext_eq<test_type>(*it, *rng_it);
}

TYPED_TEST_P(iterator_fixture, move_forward_pre)
{
    move_forward_pre_test<TypeParam>(this->begin_iterator, this->end_iterator, this->range_to_compare);

    if constexpr (!std::Same<typename TestFixture::iterator_tag, std::input_iterator_tag>) // iterate over it again
        move_forward_pre_test<TypeParam>(this->begin_iterator, this->end_iterator, this->range_to_compare);

    if constexpr (TestFixture::const_iterable)
        move_forward_pre_test<TypeParam>(this->begin_const_iterator, this->end_const_iterator, this->range_to_compare);
}

TYPED_TEST_P(iterator_fixture, move_forward_post)
{
    move_forward_post_test<TypeParam>(this->begin_iterator, this->end_iterator, this->range_to_compare);

    if constexpr (!std::Same<typename TestFixture::iterator_tag, std::input_iterator_tag>) // iterate over it again
        move_forward_post_test<TypeParam>(this->begin_iterator, this->end_iterator, this->range_to_compare);

    if constexpr (TestFixture::const_iterable)
        move_forward_post_test<TypeParam>(this->begin_const_iterator, this->end_const_iterator, this->range_to_compare);
}

// ---------------------------------------------------------------------------------------------------------------------
// Bidirectional Iterator
// ---------------------------------------------------------------------------------------------------------------------

template <typename test_type, typename it_begin_t, typename it_sentinel_t, typename rng_t>
inline void move_backward_test(it_begin_t && it_begin, it_sentinel_t && it_end, rng_t && rng)
{
    // move to last position
    auto pre_end_it = it_begin;
    auto rng_pre_end_it = rng.begin();

    for (; std::ranges::next(pre_end_it) != it_end; ++pre_end_it, ++rng_pre_end_it);

    // pre-decrement
    auto rng_it = rng_pre_end_it;
    for (auto it = pre_end_it; it != it_begin; --it, --rng_it)
        expext_eq<test_type>(*it, *rng_it);

    // post-decrement
    rng_it = rng_pre_end_it;
    for (auto it = pre_end_it; it != it_begin; --rng_it)
        expext_eq<test_type>(*(it--), *rng_it);
}

TYPED_TEST_P(iterator_fixture, move_backward)
{
    if constexpr (std::Same<typename TestFixture::iterator_tag, std::bidirectional_iterator_tag>)
    {
        move_backward_test<TypeParam>(this->begin_iterator, this->end_iterator, this->range_to_compare);

        if constexpr (TestFixture::const_iterable)
            move_backward_test<TypeParam>(this->begin_const_iterator, this->end_const_iterator, this->range_to_compare);
    }
}

// ---------------------------------------------------------------------------------------------------------------------
// Random Access Iterator
// ---------------------------------------------------------------------------------------------------------------------

template <typename test_type, typename it_begin_t, typename rng_t>
inline void jump_forward_test(it_begin_t && it_begin, rng_t && rng)
{
    size_t sz = std::ranges::distance(rng);

    // Forward
    for (size_t n = 0; n < sz; ++n)
    {
        auto it = it_begin;
        expext_eq<test_type>(rng[n], *(it += n));
        expext_eq<test_type>(rng[n], *(it));
    }

    // Forward copy
    for (size_t n = 0; n < sz; ++n)
    {
        expext_eq<test_type>(rng[n], *(it_begin + n));
        expext_eq<test_type>(rng[0], *it_begin);
    }

    // Forward copy friend
    for (size_t n = 0; n < sz; ++n)
    {
        expext_eq<test_type>(rng[n], *(n + it_begin));
        expext_eq<test_type>(rng[0], *it_begin);
    }
}

TYPED_TEST_P(iterator_fixture, jump_forward)
{
    if constexpr (std::Same<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        jump_forward_test<TypeParam>(this->begin_iterator, this->range_to_compare);

        if constexpr (TestFixture::const_iterable)
            jump_forward_test<TypeParam>(this->begin_const_iterator, this->range_to_compare);
    }
}

template <typename test_type, typename it_begin_t, typename rng_t>
inline void jump_backward_test(it_begin_t && it_begin, rng_t && rng)
{
    size_t sz = std::ranges::distance(rng);

    auto pre_end_it = it_begin + sz - 1;

    // Backward
    for (size_t n = 0; n < sz; ++n)
    {
        auto it = pre_end_it;
        expext_eq<test_type>(rng[sz - 1 - n], *(it -= n));
        expext_eq<test_type>(rng[sz - 1 - n], *it);
    }

    // Backward copy
    for (size_t n = 0; n < sz; ++n)
    {
        expext_eq<test_type>(rng[sz - n - 1], *(pre_end_it - n));
        expext_eq<test_type>(rng[sz - 1], *pre_end_it);
    }

    // Backward copy friend through (-n) + it
    for (size_t n = 0; n < sz; ++n)
    {
        expext_eq<test_type>(rng[sz - n - 1], *((-1 * n) + pre_end_it));
        expext_eq<test_type>(rng[sz - 1], *pre_end_it);
    }
}

TYPED_TEST_P(iterator_fixture, jump_backward)
{
    if constexpr (std::Same<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        jump_backward_test<TypeParam>(this->begin_iterator, this->range_to_compare);

        if constexpr (TestFixture::const_iterable)
            jump_backward_test<TypeParam>(this->begin_const_iterator, this->range_to_compare);
    }
}

template <typename test_type, typename it_begin_t, typename rng_t>
inline void jump_random_test(it_begin_t && it_begin, rng_t && rng)
{
    size_t sz = std::ranges::distance(rng);

    for (size_t n = 0; n < sz; ++n)
        expext_eq<test_type>(rng[n], it_begin[n]);
}

TYPED_TEST_P(iterator_fixture, jump_random)
{
    if constexpr (std::Same<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        jump_random_test<TypeParam>(this->begin_iterator, this->range_to_compare);

        if constexpr (TestFixture::const_iterable)
            jump_random_test<TypeParam>(this->begin_const_iterator, this->range_to_compare);
    }
}

template <typename it_begin_t, typename rng_t>
inline void difference_test(it_begin_t && it_begin, rng_t && rng)
{
    size_t sz = std::ranges::distance(rng);
    using difference_t = typename std::iterator_traits<std::remove_reference_t<decltype(it_begin)>>::difference_type;

    for (size_t n = 0; n < sz; ++n)
        EXPECT_EQ(static_cast<difference_t>(n), ((it_begin + n) - it_begin));
}

TYPED_TEST_P(iterator_fixture, difference)
{
    if constexpr (std::Same<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        difference_test(this->begin_iterator, this->range_to_compare);

        if constexpr (TestFixture::const_iterable)
            difference_test(this->begin_const_iterator, this->range_to_compare);
    }
}

TYPED_TEST_P(iterator_fixture, compare_less)
{
    if constexpr (std::Same<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        EXPECT_FALSE(this->begin_iterator < this->begin_iterator);
        EXPECT_TRUE(this->begin_iterator < this->end_iterator);
        EXPECT_TRUE(this->begin_iterator < std::ranges::next(this->begin_iterator));
    }

    if constexpr (std::Same<typename TestFixture::iterator_tag, std::random_access_iterator_tag> &&
                  TestFixture::const_iterable)
    {
        EXPECT_FALSE(this->begin_const_iterator < this->begin_const_iterator);
        EXPECT_TRUE(this->begin_const_iterator < this->end_const_iterator);
        EXPECT_TRUE(this->begin_const_iterator < std::ranges::next(this->begin_const_iterator));

        // mix
        EXPECT_FALSE(this->begin_iterator < this->begin_const_iterator);
        EXPECT_TRUE(this->begin_iterator < this->end_const_iterator);
        EXPECT_TRUE(this->begin_iterator < std::ranges::next(this->begin_const_iterator));
        EXPECT_FALSE(this->begin_const_iterator < this->begin_iterator);
        EXPECT_TRUE(this->begin_const_iterator < this->end_iterator);
        EXPECT_TRUE(this->begin_const_iterator < std::ranges::next(this->begin_iterator));
    }
}

TYPED_TEST_P(iterator_fixture, compare_greater)
{
    if constexpr (std::Same<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        EXPECT_FALSE(this->begin_iterator > this->begin_iterator);
        EXPECT_TRUE(this->end_iterator > this->begin_iterator);
        EXPECT_FALSE(this->begin_iterator > std::ranges::next(this->begin_iterator));
    }

    if constexpr (std::Same<typename TestFixture::iterator_tag, std::random_access_iterator_tag> &&
                  TestFixture::const_iterable)
    {
        EXPECT_FALSE(this->begin_const_iterator > this->begin_const_iterator);
        EXPECT_TRUE(this->end_const_iterator > this->begin_const_iterator);
        EXPECT_FALSE(this->begin_const_iterator > std::ranges::next(this->begin_const_iterator));

        // mix
        EXPECT_FALSE(this->begin_iterator > this->begin_const_iterator);
        EXPECT_TRUE(this->end_iterator > this->begin_const_iterator);
        EXPECT_FALSE(this->begin_iterator > std::ranges::next(this->begin_const_iterator));
        EXPECT_FALSE(this->begin_const_iterator > this->begin_iterator);
        EXPECT_TRUE(this->end_const_iterator > this->begin_iterator);
        EXPECT_FALSE(this->begin_const_iterator > std::ranges::next(this->begin_iterator));
    }
}

TYPED_TEST_P(iterator_fixture, compare_leq)
{
    if constexpr (std::Same<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        EXPECT_TRUE(this->begin_iterator <= this->begin_iterator);
        EXPECT_TRUE(this->begin_iterator <= this->end_iterator);
        EXPECT_TRUE(this->begin_iterator <= std::ranges::next(this->begin_iterator));
    }

    if constexpr (std::Same<typename TestFixture::iterator_tag, std::random_access_iterator_tag> &&
                  TestFixture::const_iterable)
    {
        EXPECT_TRUE(this->begin_const_iterator <= this->begin_const_iterator);
        EXPECT_TRUE(this->begin_const_iterator <= this->end_const_iterator);
        EXPECT_TRUE(this->begin_const_iterator <= std::ranges::next(this->begin_const_iterator));

        // mix
        EXPECT_TRUE(this->begin_iterator <= this->begin_const_iterator);
        EXPECT_TRUE(this->begin_iterator <= this->end_const_iterator);
        EXPECT_TRUE(this->begin_iterator <= std::ranges::next(this->begin_const_iterator));
        EXPECT_TRUE(this->begin_const_iterator <= this->begin_iterator);
        EXPECT_TRUE(this->begin_const_iterator <= this->end_iterator);
        EXPECT_TRUE(this->begin_const_iterator <= std::ranges::next(this->begin_iterator));
    }
}

TYPED_TEST_P(iterator_fixture, compare_geq)
{
    if constexpr (std::Same<typename TestFixture::iterator_tag, std::random_access_iterator_tag>)
    {
        EXPECT_TRUE(this->begin_iterator >= this->begin_iterator);
        EXPECT_TRUE(this->end_iterator >= this->begin_iterator);
        EXPECT_FALSE(this->begin_iterator >= std::ranges::next(this->begin_iterator));
    }

    if constexpr (std::Same<typename TestFixture::iterator_tag, std::random_access_iterator_tag> &&
                  TestFixture::const_iterable)
    {
        EXPECT_TRUE(this->begin_const_iterator >= this->begin_const_iterator);
        EXPECT_TRUE(this->end_const_iterator >= this->begin_const_iterator);
        EXPECT_FALSE(this->begin_const_iterator >= std::ranges::next(this->begin_const_iterator));

        // mix
        EXPECT_TRUE(this->begin_iterator >= this->begin_const_iterator);
        EXPECT_TRUE(this->end_iterator >= this->begin_const_iterator);
        EXPECT_FALSE(this->begin_iterator >= std::ranges::next(this->begin_const_iterator));
        EXPECT_TRUE(this->begin_const_iterator >= this->begin_iterator);
        EXPECT_TRUE(this->end_const_iterator >= this->begin_iterator);
        EXPECT_FALSE(this->begin_const_iterator >= std::ranges::next(this->begin_iterator));
    }
}

REGISTER_TYPED_TEST_CASE_P(iterator_fixture,
                           concept_check,
                           const_non_const_compatibility,
                           dereference,
                           compare,
                           move_forward_pre,
                           move_forward_post,
                           move_backward,
                           jump_forward,
                           jump_backward,
                           jump_random,
                           difference,
                           compare_less,
                           compare_greater,
                           compare_leq,
                           compare_geq);
