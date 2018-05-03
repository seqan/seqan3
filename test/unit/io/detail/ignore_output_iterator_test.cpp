// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
// ==========================================================================

#include <gtest/gtest.h>

#include <seqan3/std/concept/iterator.hpp>
#include <seqan3/io/detail/ignore_output_iterator.hpp>

using namespace seqan3;
using namespace std::literals;

TEST(ignore_output_iterator, concept)
{
    EXPECT_TRUE((output_iterator_concept<detail::ignore_output_iterator, char>));
    EXPECT_TRUE((output_iterator_concept<detail::ignore_output_iterator, int>));
    EXPECT_FALSE((input_iterator_concept<detail::ignore_output_iterator>));
}

TEST(ignore_output_iterator, assign)
{
    detail::ignore_output_iterator it;
    EXPECT_EQ(std::addressof(it = 'A'), std::addressof(it));
    EXPECT_EQ(std::addressof(it = 10), std::addressof(it));
}

TEST(ignore_output_iterator, pre_increment)
{
    detail::ignore_output_iterator it;
    EXPECT_EQ(std::addressof(++it), std::addressof(it));
}

TEST(ignore_output_iterator, post_increment)
{
    detail::ignore_output_iterator it;
    EXPECT_EQ(std::addressof(it++), std::addressof(it));
}

TEST(ignore_output_iterator, dereference)
{
    detail::ignore_output_iterator it;
    EXPECT_EQ(std::addressof(*it), std::addressof(it));
}

TEST(ignore_output_iterator, make_conversion_output_iterator)
{
    auto it = detail::make_conversion_output_iterator(std::ignore);
    EXPECT_TRUE((std::is_same_v<decltype(it), detail::ignore_output_iterator>));
}
