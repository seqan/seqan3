// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include "auxiliary_iterator.hpp"

TEST(iterator_concepts, readable)
{
    EXPECT_TRUE((std::readable<input_iterator>));
    EXPECT_TRUE((!std::readable<output_iterator>));
    EXPECT_TRUE((std::readable<forward_iterator>));
    EXPECT_TRUE((std::readable<bidirectional_iterator>));
    EXPECT_TRUE((std::readable<random_access_iterator>));
    EXPECT_TRUE((std::readable<forward_iterator_const>));
    EXPECT_TRUE((std::readable<bidirectional_iterator_const>));
    EXPECT_TRUE((std::readable<random_access_iterator_const>));
}

TEST(iterator_concepts, writable)
{
    EXPECT_TRUE((!std::writable<input_iterator, char>));
    EXPECT_TRUE((std::writable<output_iterator, char>));
    EXPECT_TRUE((std::writable<forward_iterator, char>));
    EXPECT_TRUE((!std::writable<forward_iterator_const, char>));
    EXPECT_TRUE((std::writable<bidirectional_iterator, char>));
    EXPECT_TRUE((!std::writable<bidirectional_iterator_const, char>));
    EXPECT_TRUE((std::writable<random_access_iterator, char>));
    EXPECT_TRUE((!std::writable<random_access_iterator_const, char>));
}

TEST(iterator_concepts, weakly_incrementable)
{
    EXPECT_TRUE((std::weakly_incrementable<input_iterator>));
    EXPECT_TRUE((std::weakly_incrementable<output_iterator>));
    EXPECT_TRUE((std::weakly_incrementable<forward_iterator>));
    EXPECT_TRUE((std::weakly_incrementable<forward_iterator_const>));
    EXPECT_TRUE((std::weakly_incrementable<bidirectional_iterator>));
    EXPECT_TRUE((std::weakly_incrementable<bidirectional_iterator_const>));
    EXPECT_TRUE((std::weakly_incrementable<random_access_iterator>));
    EXPECT_TRUE((std::weakly_incrementable<random_access_iterator_const>));
}

TEST(iterator_concepts, incrementable)
{
    EXPECT_TRUE((std::incrementable<input_iterator>));
    EXPECT_TRUE((!std::incrementable<output_iterator>));
    EXPECT_TRUE((std::incrementable<forward_iterator>));
    EXPECT_TRUE((std::incrementable<forward_iterator_const>));
    EXPECT_TRUE((std::incrementable<bidirectional_iterator>));
    EXPECT_TRUE((std::incrementable<bidirectional_iterator_const>));
    EXPECT_TRUE((std::incrementable<random_access_iterator>));
    EXPECT_TRUE((std::incrementable<random_access_iterator_const>));
}

TEST(iterator_concepts, Iterator)
{
    EXPECT_TRUE((std::input_or_output_iterator<input_iterator>));
    EXPECT_TRUE((std::input_or_output_iterator<output_iterator>));
    EXPECT_TRUE((std::input_or_output_iterator<forward_iterator>));
    EXPECT_TRUE((std::input_or_output_iterator<forward_iterator_const>));
    EXPECT_TRUE((std::input_or_output_iterator<bidirectional_iterator>));
    EXPECT_TRUE((std::input_or_output_iterator<bidirectional_iterator_const>));
    EXPECT_TRUE((std::input_or_output_iterator<random_access_iterator>));
    EXPECT_TRUE((std::input_or_output_iterator<random_access_iterator_const>));
}

TEST(iterator_concepts, sentinel_for)
{
    EXPECT_TRUE((std::sentinel_for<test_sentinel<char>,
                                          input_iterator>));
    EXPECT_TRUE((std::sentinel_for<test_sentinel<char>,
                                          output_iterator>));
    EXPECT_TRUE((std::sentinel_for<test_sentinel<char>,
                                          forward_iterator>));
    EXPECT_TRUE((std::sentinel_for<test_sentinel<char>,
                                          forward_iterator_const>));
    EXPECT_TRUE((std::sentinel_for<test_sentinel<char>,
                                          bidirectional_iterator>));
    EXPECT_TRUE((std::sentinel_for<test_sentinel<char>,
                                          bidirectional_iterator_const>));
    EXPECT_TRUE((std::sentinel_for<test_sentinel<char>,
                                          random_access_iterator>));
    EXPECT_TRUE((std::sentinel_for<test_sentinel<char>,
                                          random_access_iterator_const>));

    EXPECT_TRUE((std::sentinel_for<test_sized_sentinel<
                                          input_iterator>,
                                          input_iterator>));
    EXPECT_TRUE((std::is_same_v<char, value_type_t<seqan3::ostream_iterator<char>>>));
    EXPECT_TRUE((std::sentinel_for<test_sized_sentinel<
                                          output_iterator>,
                                          output_iterator>));
    EXPECT_TRUE((std::sentinel_for<test_sized_sentinel<
                                          forward_iterator>,
                                          forward_iterator>));
    EXPECT_TRUE((std::sentinel_for<test_sized_sentinel<
                                          forward_iterator_const>,
                                          forward_iterator_const>));
    EXPECT_TRUE((std::sentinel_for<test_sized_sentinel<
                                          bidirectional_iterator>,
                                          bidirectional_iterator>));
    EXPECT_TRUE((std::sentinel_for<test_sized_sentinel<
                                          bidirectional_iterator_const>,
                                          bidirectional_iterator_const>));
    EXPECT_TRUE((std::sentinel_for<test_sized_sentinel<
                                          random_access_iterator>,
                                          random_access_iterator>));
    EXPECT_TRUE((std::sentinel_for<test_sized_sentinel<
                                          random_access_iterator_const>,
                                          random_access_iterator_const>));
}

TEST(iterator_concepts, sized_sentinel_for)
{
    EXPECT_TRUE((!std::sized_sentinel_for<test_sentinel<char>,
                                                 input_iterator>));
    EXPECT_TRUE((!std::sized_sentinel_for<test_sentinel<char>,
                                                 output_iterator>));
    EXPECT_TRUE((!std::sized_sentinel_for<test_sentinel<char>,
                                                 forward_iterator>));
    EXPECT_TRUE((!std::sized_sentinel_for<test_sentinel<char>,
                                                 forward_iterator_const>));
    EXPECT_TRUE((!std::sized_sentinel_for<test_sentinel<char>,
                                                 bidirectional_iterator>));
    EXPECT_TRUE((!std::sized_sentinel_for<test_sentinel<char>,
                                                 bidirectional_iterator_const>));
    EXPECT_TRUE((!std::sized_sentinel_for<test_sentinel<char>,
                                                 random_access_iterator>));
    EXPECT_TRUE((!std::sized_sentinel_for<test_sentinel<char>,
                                                 random_access_iterator_const>));

    EXPECT_TRUE((!std::sized_sentinel_for<test_sized_sentinel<
                                                 input_iterator>,
                                                 input_iterator>));
    EXPECT_TRUE((!std::sized_sentinel_for<test_sized_sentinel<
                                                 output_iterator>,
                                                 output_iterator>));
    EXPECT_TRUE((!std::sized_sentinel_for<test_sized_sentinel<
                                                 forward_iterator>,
                                                 forward_iterator>));
    EXPECT_TRUE((!std::sized_sentinel_for<test_sized_sentinel<
                                                 forward_iterator_const>,
                                                 forward_iterator_const>));
    EXPECT_TRUE((!std::sized_sentinel_for<test_sized_sentinel<
                                                 bidirectional_iterator>,
                                                 bidirectional_iterator>));
    EXPECT_TRUE((!std::sized_sentinel_for<test_sized_sentinel<
                                                 bidirectional_iterator_const>,
                                                 bidirectional_iterator_const>));
    EXPECT_TRUE((std::sized_sentinel_for<test_sized_sentinel<
                                                random_access_iterator>,
                                                random_access_iterator>));
    EXPECT_TRUE((std::sized_sentinel_for<test_sized_sentinel<
                                                random_access_iterator_const>,
                                                random_access_iterator_const>));
}

TEST(iterator_concepts, output_iterator)
{
    EXPECT_TRUE((!std::output_iterator<input_iterator, char>));
    EXPECT_TRUE((std::output_iterator<output_iterator, char>));
    EXPECT_TRUE((std::output_iterator<forward_iterator, char>));
    EXPECT_TRUE((!std::output_iterator<forward_iterator_const, char>));
    EXPECT_TRUE((std::output_iterator<bidirectional_iterator, char>));
    EXPECT_TRUE((!std::output_iterator<bidirectional_iterator_const, char>));
    EXPECT_TRUE((std::output_iterator<random_access_iterator, char>));
    EXPECT_TRUE((!std::output_iterator<random_access_iterator_const, char>));
}

TEST(iterator_concepts, input_iterator)
{
    EXPECT_TRUE((std::input_iterator<input_iterator>));
    EXPECT_TRUE((!std::input_iterator<output_iterator>));
    EXPECT_TRUE((std::input_iterator<forward_iterator>));
    EXPECT_TRUE((std::input_iterator<forward_iterator_const>));
    EXPECT_TRUE((std::input_iterator<bidirectional_iterator>));
    EXPECT_TRUE((std::input_iterator<bidirectional_iterator_const>));
    EXPECT_TRUE((std::input_iterator<random_access_iterator>));
    EXPECT_TRUE((std::input_iterator<random_access_iterator_const>));
}

TEST(iterator_concepts, forward_iterator)
{
    EXPECT_TRUE((!std::forward_iterator<input_iterator>));
    EXPECT_TRUE((!std::forward_iterator<output_iterator>));
    EXPECT_TRUE((std::forward_iterator<forward_iterator>));
    EXPECT_TRUE((std::forward_iterator<forward_iterator_const>));
    EXPECT_TRUE((std::forward_iterator<bidirectional_iterator>));
    EXPECT_TRUE((std::forward_iterator<bidirectional_iterator_const>));
    EXPECT_TRUE((std::forward_iterator<random_access_iterator>));
    EXPECT_TRUE((std::forward_iterator<random_access_iterator_const>));
}

TEST(iterator_concepts, bidirectional_iterator)
{
    EXPECT_TRUE((!std::bidirectional_iterator<input_iterator>));
    EXPECT_TRUE((!std::bidirectional_iterator<output_iterator>));
    EXPECT_TRUE((!std::bidirectional_iterator<forward_iterator>));
    EXPECT_TRUE((!std::bidirectional_iterator<forward_iterator_const>));
    EXPECT_TRUE((std::bidirectional_iterator<bidirectional_iterator>));
    EXPECT_TRUE((std::bidirectional_iterator<bidirectional_iterator_const>));
    EXPECT_TRUE((std::bidirectional_iterator<random_access_iterator>));
    EXPECT_TRUE((std::bidirectional_iterator<random_access_iterator_const>));
}

TEST(iterator_concepts, random_access_iterator)
{
    EXPECT_TRUE((!std::random_access_iterator<input_iterator>));
    EXPECT_TRUE((!std::random_access_iterator<output_iterator>));
    EXPECT_TRUE((!std::random_access_iterator<forward_iterator>));
    EXPECT_TRUE((!std::random_access_iterator<forward_iterator_const>));
    EXPECT_TRUE((!std::random_access_iterator<bidirectional_iterator>));
    EXPECT_TRUE((!std::random_access_iterator<bidirectional_iterator_const>));
    EXPECT_TRUE((std::random_access_iterator<random_access_iterator>));
    EXPECT_TRUE((std::random_access_iterator<random_access_iterator_const>));
}
