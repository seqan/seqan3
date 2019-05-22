// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include "auxiliary_iterator.hpp"

TEST(iterator_concepts, Readable)
{
    EXPECT_TRUE((std::Readable<input_iterator>));
    EXPECT_TRUE((!std::Readable<output_iterator>));
    EXPECT_TRUE((std::Readable<forward_iterator>));
    EXPECT_TRUE((std::Readable<bidirectional_iterator>));
    EXPECT_TRUE((std::Readable<random_access_iterator>));
    EXPECT_TRUE((std::Readable<forward_iterator_const>));
    EXPECT_TRUE((std::Readable<bidirectional_iterator_const>));
    EXPECT_TRUE((std::Readable<random_access_iterator_const>));
}

TEST(iterator_concepts, Writable)
{
    EXPECT_TRUE((!std::Writable<input_iterator, char>));
    EXPECT_TRUE((std::Writable<output_iterator, char>));
    EXPECT_TRUE((std::Writable<forward_iterator, char>));
    EXPECT_TRUE((!std::Writable<forward_iterator_const, char>));
    EXPECT_TRUE((std::Writable<bidirectional_iterator, char>));
    EXPECT_TRUE((!std::Writable<bidirectional_iterator_const, char>));
    EXPECT_TRUE((std::Writable<random_access_iterator, char>));
    EXPECT_TRUE((!std::Writable<random_access_iterator_const, char>));
}

TEST(iterator_concepts, WeaklyIncrementable)
{
    EXPECT_TRUE((std::WeaklyIncrementable<input_iterator>));
    EXPECT_TRUE((std::WeaklyIncrementable<output_iterator>));
    EXPECT_TRUE((std::WeaklyIncrementable<forward_iterator>));
    EXPECT_TRUE((std::WeaklyIncrementable<forward_iterator_const>));
    EXPECT_TRUE((std::WeaklyIncrementable<bidirectional_iterator>));
    EXPECT_TRUE((std::WeaklyIncrementable<bidirectional_iterator_const>));
    EXPECT_TRUE((std::WeaklyIncrementable<random_access_iterator>));
    EXPECT_TRUE((std::WeaklyIncrementable<random_access_iterator_const>));
}

TEST(iterator_concepts, Incrementable)
{
    EXPECT_TRUE((std::Incrementable<input_iterator>));
    EXPECT_TRUE((!std::Incrementable<output_iterator>));
    EXPECT_TRUE((std::Incrementable<forward_iterator>));
    EXPECT_TRUE((std::Incrementable<forward_iterator_const>));
    EXPECT_TRUE((std::Incrementable<bidirectional_iterator>));
    EXPECT_TRUE((std::Incrementable<bidirectional_iterator_const>));
    EXPECT_TRUE((std::Incrementable<random_access_iterator>));
    EXPECT_TRUE((std::Incrementable<random_access_iterator_const>));
}

TEST(iterator_concepts, Iterator)
{
    EXPECT_TRUE((std::Iterator<input_iterator>));
    EXPECT_TRUE((std::Iterator<output_iterator>));
    EXPECT_TRUE((std::Iterator<forward_iterator>));
    EXPECT_TRUE((std::Iterator<forward_iterator_const>));
    EXPECT_TRUE((std::Iterator<bidirectional_iterator>));
    EXPECT_TRUE((std::Iterator<bidirectional_iterator_const>));
    EXPECT_TRUE((std::Iterator<random_access_iterator>));
    EXPECT_TRUE((std::Iterator<random_access_iterator_const>));
}

TEST(iterator_concepts, Sentinel)
{
    EXPECT_TRUE((std::Sentinel<test_sentinel<char>,
                                          input_iterator>));
    EXPECT_TRUE((std::Sentinel<test_sentinel<char>,
                                          output_iterator>));
    EXPECT_TRUE((std::Sentinel<test_sentinel<char>,
                                          forward_iterator>));
    EXPECT_TRUE((std::Sentinel<test_sentinel<char>,
                                          forward_iterator_const>));
    EXPECT_TRUE((std::Sentinel<test_sentinel<char>,
                                          bidirectional_iterator>));
    EXPECT_TRUE((std::Sentinel<test_sentinel<char>,
                                          bidirectional_iterator_const>));
    EXPECT_TRUE((std::Sentinel<test_sentinel<char>,
                                          random_access_iterator>));
    EXPECT_TRUE((std::Sentinel<test_sentinel<char>,
                                          random_access_iterator_const>));

    EXPECT_TRUE((std::Sentinel<test_sized_sentinel<
                                          input_iterator>,
                                          input_iterator>));
    EXPECT_TRUE((std::is_same_v<char, value_type_t<seqan3::ostream_iterator<char>>>));
    EXPECT_TRUE((std::Sentinel<test_sized_sentinel<
                                          output_iterator>,
                                          output_iterator>));
    EXPECT_TRUE((std::Sentinel<test_sized_sentinel<
                                          forward_iterator>,
                                          forward_iterator>));
    EXPECT_TRUE((std::Sentinel<test_sized_sentinel<
                                          forward_iterator_const>,
                                          forward_iterator_const>));
    EXPECT_TRUE((std::Sentinel<test_sized_sentinel<
                                          bidirectional_iterator>,
                                          bidirectional_iterator>));
    EXPECT_TRUE((std::Sentinel<test_sized_sentinel<
                                          bidirectional_iterator_const>,
                                          bidirectional_iterator_const>));
    EXPECT_TRUE((std::Sentinel<test_sized_sentinel<
                                          random_access_iterator>,
                                          random_access_iterator>));
    EXPECT_TRUE((std::Sentinel<test_sized_sentinel<
                                          random_access_iterator_const>,
                                          random_access_iterator_const>));
}

TEST(iterator_concepts, SizedSentinel)
{
    EXPECT_TRUE((!std::SizedSentinel<test_sentinel<char>,
                                                 input_iterator>));
    EXPECT_TRUE((!std::SizedSentinel<test_sentinel<char>,
                                                 output_iterator>));
    EXPECT_TRUE((!std::SizedSentinel<test_sentinel<char>,
                                                 forward_iterator>));
    EXPECT_TRUE((!std::SizedSentinel<test_sentinel<char>,
                                                 forward_iterator_const>));
    EXPECT_TRUE((!std::SizedSentinel<test_sentinel<char>,
                                                 bidirectional_iterator>));
    EXPECT_TRUE((!std::SizedSentinel<test_sentinel<char>,
                                                 bidirectional_iterator_const>));
    EXPECT_TRUE((!std::SizedSentinel<test_sentinel<char>,
                                                 random_access_iterator>));
    EXPECT_TRUE((!std::SizedSentinel<test_sentinel<char>,
                                                 random_access_iterator_const>));

    EXPECT_TRUE((!std::SizedSentinel<test_sized_sentinel<
                                                 input_iterator>,
                                                 input_iterator>));
    EXPECT_TRUE((!std::SizedSentinel<test_sized_sentinel<
                                                 output_iterator>,
                                                 output_iterator>));
    EXPECT_TRUE((!std::SizedSentinel<test_sized_sentinel<
                                                 forward_iterator>,
                                                 forward_iterator>));
    EXPECT_TRUE((!std::SizedSentinel<test_sized_sentinel<
                                                 forward_iterator_const>,
                                                 forward_iterator_const>));
    EXPECT_TRUE((!std::SizedSentinel<test_sized_sentinel<
                                                 bidirectional_iterator>,
                                                 bidirectional_iterator>));
    EXPECT_TRUE((!std::SizedSentinel<test_sized_sentinel<
                                                 bidirectional_iterator_const>,
                                                 bidirectional_iterator_const>));
    EXPECT_TRUE((std::SizedSentinel<test_sized_sentinel<
                                                random_access_iterator>,
                                                random_access_iterator>));
    EXPECT_TRUE((std::SizedSentinel<test_sized_sentinel<
                                                random_access_iterator_const>,
                                                random_access_iterator_const>));
}

TEST(iterator_concepts, OutputIterator)
{
    EXPECT_TRUE((!std::OutputIterator<input_iterator, char>));
    EXPECT_TRUE((std::OutputIterator<output_iterator, char>));
    EXPECT_TRUE((std::OutputIterator<forward_iterator, char>));
    EXPECT_TRUE((!std::OutputIterator<forward_iterator_const, char>));
    EXPECT_TRUE((std::OutputIterator<bidirectional_iterator, char>));
    EXPECT_TRUE((!std::OutputIterator<bidirectional_iterator_const, char>));
    EXPECT_TRUE((std::OutputIterator<random_access_iterator, char>));
    EXPECT_TRUE((!std::OutputIterator<random_access_iterator_const, char>));
}

TEST(iterator_concepts, InputIterator)
{
    EXPECT_TRUE((std::InputIterator<input_iterator>));
    EXPECT_TRUE((!std::InputIterator<output_iterator>));
    EXPECT_TRUE((std::InputIterator<forward_iterator>));
    EXPECT_TRUE((std::InputIterator<forward_iterator_const>));
    EXPECT_TRUE((std::InputIterator<bidirectional_iterator>));
    EXPECT_TRUE((std::InputIterator<bidirectional_iterator_const>));
    EXPECT_TRUE((std::InputIterator<random_access_iterator>));
    EXPECT_TRUE((std::InputIterator<random_access_iterator_const>));
}

TEST(iterator_concepts, ForwardIterator)
{
    EXPECT_TRUE((!std::ForwardIterator<input_iterator>));
    EXPECT_TRUE((!std::ForwardIterator<output_iterator>));
    EXPECT_TRUE((std::ForwardIterator<forward_iterator>));
    EXPECT_TRUE((std::ForwardIterator<forward_iterator_const>));
    EXPECT_TRUE((std::ForwardIterator<bidirectional_iterator>));
    EXPECT_TRUE((std::ForwardIterator<bidirectional_iterator_const>));
    EXPECT_TRUE((std::ForwardIterator<random_access_iterator>));
    EXPECT_TRUE((std::ForwardIterator<random_access_iterator_const>));
}

TEST(iterator_concepts, BidirectionalIterator)
{
    EXPECT_TRUE((!std::BidirectionalIterator<input_iterator>));
    EXPECT_TRUE((!std::BidirectionalIterator<output_iterator>));
    EXPECT_TRUE((!std::BidirectionalIterator<forward_iterator>));
    EXPECT_TRUE((!std::BidirectionalIterator<forward_iterator_const>));
    EXPECT_TRUE((std::BidirectionalIterator<bidirectional_iterator>));
    EXPECT_TRUE((std::BidirectionalIterator<bidirectional_iterator_const>));
    EXPECT_TRUE((std::BidirectionalIterator<random_access_iterator>));
    EXPECT_TRUE((std::BidirectionalIterator<random_access_iterator_const>));
}

TEST(iterator_concepts, RandomAccessIterator)
{
    EXPECT_TRUE((!std::RandomAccessIterator<input_iterator>));
    EXPECT_TRUE((!std::RandomAccessIterator<output_iterator>));
    EXPECT_TRUE((!std::RandomAccessIterator<forward_iterator>));
    EXPECT_TRUE((!std::RandomAccessIterator<forward_iterator_const>));
    EXPECT_TRUE((!std::RandomAccessIterator<bidirectional_iterator>));
    EXPECT_TRUE((!std::RandomAccessIterator<bidirectional_iterator_const>));
    EXPECT_TRUE((std::RandomAccessIterator<random_access_iterator>));
    EXPECT_TRUE((std::RandomAccessIterator<random_access_iterator_const>));
}
