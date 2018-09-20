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

<<<<<<< HEAD
#include <seqan3/io/detail/in_file_iterator.hpp>
#include <seqan3/std/concept/iterator.hpp>
=======
#include <seqan3/std/iterator>

#include <seqan3/io/detail/in_file_iterator.hpp>
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6

using namespace seqan3;

//NOTE(h-2): This class is extensively tested via *_file_in. This is just a minimal test.

struct fake_file_t : std::vector<int>
{
    using base = std::vector<int>;

    using base::base;

    size_t current_position = 0;
    bool at_end = false;

    void read_next_record()
    {
        ++current_position;
        if (current_position == size())
            at_end = true;
    }

    int & front()
    {
        return begin()[current_position];
    }

    fake_file_t() :
        base{}
    {}
};

TEST(in_file_iterator, concepts)
{
    using it_t = detail::in_file_iterator<fake_file_t>;

<<<<<<< HEAD
    EXPECT_TRUE((input_iterator_concept<it_t>));
=======
    EXPECT_TRUE((std::InputIterator<it_t>));
>>>>>>> 41b42cc5d45c544a427ed079af957ad4366ea9e6
}

TEST(in_file_iterator, member_types)
{
    using it_t = detail::in_file_iterator<fake_file_t>;
    EXPECT_TRUE((std::is_same_v<typename it_t::value_type,
                                int>));
    EXPECT_TRUE((std::is_same_v<typename it_t::reference,
                                int &>));
    EXPECT_TRUE((std::is_same_v<typename it_t::const_reference,
                                int &>));
    EXPECT_TRUE((std::is_same_v<typename it_t::difference_type,
                                std::ptrdiff_t>));
    EXPECT_TRUE((std::is_same_v<typename it_t::size_type,
                                size_t>));
    EXPECT_TRUE((std::is_same_v<typename it_t::iterator_category,
                                std::input_iterator_tag>));
}

TEST(in_file_iterator, operations)
{
    using it_t = detail::in_file_iterator<fake_file_t>;

    fake_file_t f{1, 2, 3, 4, 5, 6, 8};

    // construct
    it_t it{f};

    // deref
    EXPECT_EQ(*it, 1);

    // pre-inc
    EXPECT_EQ(*(++it), 2);

    // post-inc
    it++;

    // deref
    EXPECT_EQ(*it, 3);
}

TEST(in_file_iterator, comparison)
{
    using it_t = detail::in_file_iterator<fake_file_t>;

    fake_file_t f{1, 2, 3, 4, 5, 6, 8};
    it_t it{f};

    // not at end
    EXPECT_FALSE(it == ranges::default_sentinel{});

    // consume the entire range
    ++it; ++it; ++it; ++it; ++it; ++it; ++it;

    // at end
    EXPECT_TRUE(it == ranges::default_sentinel{});
}
