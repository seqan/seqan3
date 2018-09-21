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

#include <vector>

#include <seqan3/std/iterator>

#include <seqan3/io/detail/out_file_iterator.hpp>

using namespace seqan3;

//NOTE(h-2): This class is extensively tested via *_file_out. This is just a minimal test.

TEST(out_file_iterator, concepts)
{
    using it_t = detail::out_file_iterator<std::vector<int>>;

    EXPECT_TRUE((std::OutputIterator<it_t, int>));
}

TEST(out_file_iterator, member_types)
{
    using it_t = detail::out_file_iterator<std::vector<int>>;
    EXPECT_TRUE((std::is_same_v<typename it_t::value_type,
                                void>));
    EXPECT_TRUE((std::is_same_v<typename it_t::reference,
                                void>));
    EXPECT_TRUE((std::is_same_v<typename it_t::const_reference,
                                void>));
    EXPECT_TRUE((std::is_same_v<typename it_t::difference_type,
                                std::ptrdiff_t>));
    EXPECT_TRUE((std::is_same_v<typename it_t::size_type,
                                void>));
    EXPECT_TRUE((std::is_same_v<typename it_t::iterator_category,
                                std::output_iterator_tag>));
}

TEST(out_file_iterator, operations)
{
    using it_t = detail::out_file_iterator<std::vector<int>>;

    std::vector<int> fake_file;

    // construct
    it_t it{fake_file};

    // pre-inc, no-op
    ++it;

    // post-inc, no-op
    it++;

    // assign to iterator
    it = 3;
    EXPECT_EQ(fake_file.size(), 1u);
    EXPECT_EQ(fake_file[0], 3);

    // assign to deref'ed iterator
    *it = 7;
    EXPECT_EQ(fake_file.size(), 2u);
    EXPECT_EQ(fake_file[0], 3);
    EXPECT_EQ(fake_file[1], 7);

    // assign to deref'ed iterator-increment
    *it++ = 9;
    EXPECT_EQ(fake_file.size(), 3u);
    EXPECT_EQ(fake_file[0], 3);
    EXPECT_EQ(fake_file[1], 7);
    EXPECT_EQ(fake_file[2], 9);
}

TEST(out_file_iterator, comparison)
{
    using it_t = detail::out_file_iterator<std::vector<int>>;

    std::vector<int> fake_file;
    it_t it{fake_file};

    // never at end
    EXPECT_FALSE(it == ranges::default_sentinel{});
}
