// ============================================================================
//                    SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//        * Redistributions of source code must retain the above copyright
//          notice, this list of conditions and the following disclaimer.
//        * Redistributions in binary form must reproduce the above copyright
//          notice, this list of conditions and the following disclaimer in the
//          documentation and/or other materials provided with the distribution.
//        * Neither the name of Knut Reinert or the FU Berlin nor the names of
//          its contributors may be used to endorse or promote products derived
//          from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

#include <gtest/gtest.h>

#include "auxiliary_iterator.hpp"

<<<<<<< HEAD
using namespace seqan3;

TEST(readable_concept, basic)
{
    EXPECT_TRUE((readable_concept<input_iterator>));
    EXPECT_TRUE((!readable_concept<output_iterator>));
    EXPECT_TRUE((readable_concept<forward_iterator>));
    EXPECT_TRUE((readable_concept<bidirectional_iterator>));
    EXPECT_TRUE((readable_concept<random_access_iterator>));
    EXPECT_TRUE((readable_concept<forward_iterator_const>));
    EXPECT_TRUE((readable_concept<bidirectional_iterator_const>));
    EXPECT_TRUE((readable_concept<random_access_iterator_const>));
}

TEST(writable_concept, basic)
{
    EXPECT_TRUE((!writable_concept<input_iterator, char>));
    EXPECT_TRUE((writable_concept<output_iterator, char>));
    EXPECT_TRUE((writable_concept<forward_iterator, char>));
    EXPECT_TRUE((!writable_concept<forward_iterator_const, char>));
    EXPECT_TRUE((writable_concept<bidirectional_iterator, char>));
    EXPECT_TRUE((!writable_concept<bidirectional_iterator_const, char>));
    EXPECT_TRUE((writable_concept<random_access_iterator, char>));
    EXPECT_TRUE((!writable_concept<random_access_iterator_const, char>));
}

TEST(weakly_incrementable_concept, basic)
{
    EXPECT_TRUE((weakly_incrementable_concept<input_iterator>));
    EXPECT_TRUE((weakly_incrementable_concept<output_iterator>));
    EXPECT_TRUE((weakly_incrementable_concept<forward_iterator>));
    EXPECT_TRUE((weakly_incrementable_concept<forward_iterator_const>));
    EXPECT_TRUE((weakly_incrementable_concept<bidirectional_iterator>));
    EXPECT_TRUE((weakly_incrementable_concept<bidirectional_iterator_const>));
    EXPECT_TRUE((weakly_incrementable_concept<random_access_iterator>));
    EXPECT_TRUE((weakly_incrementable_concept<random_access_iterator_const>));
}

TEST(incrementable_concept, basic)
{
    EXPECT_TRUE((incrementable_concept<input_iterator>));
    EXPECT_TRUE((!incrementable_concept<output_iterator>));
    EXPECT_TRUE((incrementable_concept<forward_iterator>));
    EXPECT_TRUE((incrementable_concept<forward_iterator_const>));
    EXPECT_TRUE((incrementable_concept<bidirectional_iterator>));
    EXPECT_TRUE((incrementable_concept<bidirectional_iterator_const>));
    EXPECT_TRUE((incrementable_concept<random_access_iterator>));
    EXPECT_TRUE((incrementable_concept<random_access_iterator_const>));
}

TEST(iterator_concept, basic)
{
    EXPECT_TRUE((iterator_concept<input_iterator>));
    EXPECT_TRUE((iterator_concept<output_iterator>));
    EXPECT_TRUE((iterator_concept<forward_iterator>));
    EXPECT_TRUE((iterator_concept<forward_iterator_const>));
    EXPECT_TRUE((iterator_concept<bidirectional_iterator>));
    EXPECT_TRUE((iterator_concept<bidirectional_iterator_const>));
    EXPECT_TRUE((iterator_concept<random_access_iterator>));
    EXPECT_TRUE((iterator_concept<random_access_iterator_const>));
}

TEST(sentinel_concept, basic)
{
    EXPECT_TRUE((sentinel_concept<test_sentinel<char>,
                                          input_iterator>));
    EXPECT_TRUE((sentinel_concept<test_sentinel<char>,
                                          output_iterator>));
    EXPECT_TRUE((sentinel_concept<test_sentinel<char>,
                                          forward_iterator>));
    EXPECT_TRUE((sentinel_concept<test_sentinel<char>,
                                          forward_iterator_const>));
    EXPECT_TRUE((sentinel_concept<test_sentinel<char>,
                                          bidirectional_iterator>));
    EXPECT_TRUE((sentinel_concept<test_sentinel<char>,
                                          bidirectional_iterator_const>));
    EXPECT_TRUE((sentinel_concept<test_sentinel<char>,
                                          random_access_iterator>));
    EXPECT_TRUE((sentinel_concept<test_sentinel<char>,
                                          random_access_iterator_const>));

    EXPECT_TRUE((sentinel_concept<test_sized_sentinel<
                                          input_iterator>,
                                          input_iterator>));
    EXPECT_TRUE((std::is_same_v<char, value_type_t<ranges::ostream_iterator<char>>>));
    EXPECT_TRUE((sentinel_concept<test_sized_sentinel<
                                          output_iterator>,
                                          output_iterator>));
    EXPECT_TRUE((sentinel_concept<test_sized_sentinel<
                                          forward_iterator>,
                                          forward_iterator>));
    EXPECT_TRUE((sentinel_concept<test_sized_sentinel<
                                          forward_iterator_const>,
                                          forward_iterator_const>));
    EXPECT_TRUE((sentinel_concept<test_sized_sentinel<
                                          bidirectional_iterator>,
                                          bidirectional_iterator>));
    EXPECT_TRUE((sentinel_concept<test_sized_sentinel<
                                          bidirectional_iterator_const>,
                                          bidirectional_iterator_const>));
    EXPECT_TRUE((sentinel_concept<test_sized_sentinel<
                                          random_access_iterator>,
                                          random_access_iterator>));
    EXPECT_TRUE((sentinel_concept<test_sized_sentinel<
=======
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
    EXPECT_TRUE((std::is_same_v<char, value_type_t<ranges::ostream_iterator<char>>>));
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
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6
                                          random_access_iterator_const>,
                                          random_access_iterator_const>));
}

<<<<<<< HEAD
TEST(sized_sentinel_concept, basic)
{
    EXPECT_TRUE((!sized_sentinel_concept<test_sentinel<char>,
                                                 input_iterator>));
    EXPECT_TRUE((!sized_sentinel_concept<test_sentinel<char>,
                                                 output_iterator>));
    EXPECT_TRUE((!sized_sentinel_concept<test_sentinel<char>,
                                                 forward_iterator>));
    EXPECT_TRUE((!sized_sentinel_concept<test_sentinel<char>,
                                                 forward_iterator_const>));
    EXPECT_TRUE((!sized_sentinel_concept<test_sentinel<char>,
                                                 bidirectional_iterator>));
    EXPECT_TRUE((!sized_sentinel_concept<test_sentinel<char>,
                                                 bidirectional_iterator_const>));
    EXPECT_TRUE((!sized_sentinel_concept<test_sentinel<char>,
                                                 random_access_iterator>));
    EXPECT_TRUE((!sized_sentinel_concept<test_sentinel<char>,
                                                 random_access_iterator_const>));

    EXPECT_TRUE((!sized_sentinel_concept<test_sized_sentinel<
                                                 input_iterator>,
                                                 input_iterator>));
    EXPECT_TRUE((!sized_sentinel_concept<test_sized_sentinel<
                                                 output_iterator>,
                                                 output_iterator>));
    EXPECT_TRUE((!sized_sentinel_concept<test_sized_sentinel<
                                                 forward_iterator>,
                                                 forward_iterator>));
    EXPECT_TRUE((!sized_sentinel_concept<test_sized_sentinel<
                                                 forward_iterator_const>,
                                                 forward_iterator_const>));
    EXPECT_TRUE((!sized_sentinel_concept<test_sized_sentinel<
                                                 bidirectional_iterator>,
                                                 bidirectional_iterator>));
    EXPECT_TRUE((!sized_sentinel_concept<test_sized_sentinel<
                                                 bidirectional_iterator_const>,
                                                 bidirectional_iterator_const>));
    EXPECT_TRUE((sized_sentinel_concept<test_sized_sentinel<
                                                random_access_iterator>,
                                                random_access_iterator>));
    EXPECT_TRUE((sized_sentinel_concept<test_sized_sentinel<
=======
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
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6
                                                random_access_iterator_const>,
                                                random_access_iterator_const>));
}

<<<<<<< HEAD
TEST(output_iterator_concept, basic)
{
    EXPECT_TRUE((!output_iterator_concept<input_iterator, char>));
    EXPECT_TRUE((output_iterator_concept<output_iterator, char>));
    EXPECT_TRUE((output_iterator_concept<forward_iterator, char>));
    EXPECT_TRUE((!output_iterator_concept<forward_iterator_const, char>));
    EXPECT_TRUE((output_iterator_concept<bidirectional_iterator, char>));
    EXPECT_TRUE((!output_iterator_concept<bidirectional_iterator_const, char>));
    EXPECT_TRUE((output_iterator_concept<random_access_iterator, char>));
    EXPECT_TRUE((!output_iterator_concept<random_access_iterator_const, char>));
}

TEST(input_iterator_concept, basic)
{
    EXPECT_TRUE((input_iterator_concept<input_iterator>));
    EXPECT_TRUE((!input_iterator_concept<output_iterator>));
    EXPECT_TRUE((input_iterator_concept<forward_iterator>));
    EXPECT_TRUE((input_iterator_concept<forward_iterator_const>));
    EXPECT_TRUE((input_iterator_concept<bidirectional_iterator>));
    EXPECT_TRUE((input_iterator_concept<bidirectional_iterator_const>));
    EXPECT_TRUE((input_iterator_concept<random_access_iterator>));
    EXPECT_TRUE((input_iterator_concept<random_access_iterator_const>));
}

TEST(forward_iterator_concept, basic)
{
    EXPECT_TRUE((!forward_iterator_concept<input_iterator>));
    EXPECT_TRUE((!forward_iterator_concept<output_iterator>));
    EXPECT_TRUE((forward_iterator_concept<forward_iterator>));
    EXPECT_TRUE((forward_iterator_concept<forward_iterator_const>));
    EXPECT_TRUE((forward_iterator_concept<bidirectional_iterator>));
    EXPECT_TRUE((forward_iterator_concept<bidirectional_iterator_const>));
    EXPECT_TRUE((forward_iterator_concept<random_access_iterator>));
    EXPECT_TRUE((forward_iterator_concept<random_access_iterator_const>));
}

TEST(bidirectional_iterator_concept, basic)
{
    EXPECT_TRUE((!bidirectional_iterator_concept<input_iterator>));
    EXPECT_TRUE((!bidirectional_iterator_concept<output_iterator>));
    EXPECT_TRUE((!bidirectional_iterator_concept<forward_iterator>));
    EXPECT_TRUE((!bidirectional_iterator_concept<forward_iterator_const>));
    EXPECT_TRUE((bidirectional_iterator_concept<bidirectional_iterator>));
    EXPECT_TRUE((bidirectional_iterator_concept<bidirectional_iterator_const>));
    EXPECT_TRUE((bidirectional_iterator_concept<random_access_iterator>));
    EXPECT_TRUE((bidirectional_iterator_concept<random_access_iterator_const>));
}

TEST(random_access_iterator_concept, basic)
{
    EXPECT_TRUE((!random_access_iterator_concept<input_iterator>));
    EXPECT_TRUE((!random_access_iterator_concept<output_iterator>));
    EXPECT_TRUE((!random_access_iterator_concept<forward_iterator>));
    EXPECT_TRUE((!random_access_iterator_concept<forward_iterator_const>));
    EXPECT_TRUE((!random_access_iterator_concept<bidirectional_iterator>));
    EXPECT_TRUE((!random_access_iterator_concept<bidirectional_iterator_const>));
    EXPECT_TRUE((random_access_iterator_concept<random_access_iterator>));
    EXPECT_TRUE((random_access_iterator_concept<random_access_iterator_const>));
=======
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
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6
}
