// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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

// TODO(rrahn): Make typed test for different containers?

#include <gtest/gtest.h>

#include <sstream>
#include <vector>
#include <numeric>

#include <seqan3/io/detail/direction_iterator.hpp>

using namespace seqan3;

// What we actually want is a direction_iterator over stl containers and
// And we want to make them applicable to streams
TEST(direction_iterator_test, test_direction_iterator_chunk_input_iterator_construction)
{
    using namespace seqan3::detail;

    std::vector<int> in{1, 2, 3, 4, 5, 6};

    chunk_input_iterator it{in};
    chunk_input_iterator it_2{in, true};

    EXPECT_EQ(*it, in[0]);
    EXPECT_NE(it, it_2);
}

TEST(direction_iterator_test, test_direction_iterator_chunk_input_iterator_dereference)
{
    using namespace seqan3::detail;

    std::vector<int> in{1, 2, 3, 4, 5, 6};

    chunk_input_iterator it{in};
    chunk_input_iterator const it_c{it};

    EXPECT_EQ(*it, 1);
    EXPECT_EQ(*it_c, 1);
}

TEST(direction_iterator_test, test_direction_iterator_chunk_input_iterator_increment)
{
    using namespace seqan3::detail;

    std::vector<int> in{1, 2, 3, 4, 5, 6};

    chunk_input_iterator it{in};
    chunk_input_iterator const it_c{it};

    ++it;
    EXPECT_EQ(*it++, in[1]);
    EXPECT_EQ(*it_c++, in[0]);
    ++it_c;
    EXPECT_EQ(*it_c, *it);
}

TEST(direction_iterator_test, test_direction_iterator_chunk_input_iterator_get_chunk)
{
    using namespace seqan3::detail;

    std::vector<int> in{1, 2, 3, 4, 5, 6};

    chunk_input_iterator it{in};
    chunk_input_iterator const it_c{it};

    auto rng = it.get_chunk();
    auto crng = it_c.get_chunk();
    EXPECT_EQ(*rng.begin(), in[0]);
    EXPECT_EQ(*crng.begin(), in[0]);
    EXPECT_EQ(rng.end(), std::end(in));
    EXPECT_EQ(crng.end(), std::cend(in));
}

TEST(direction_iterator_test, test_direction_iterator_chunk_input_iterator_next_chunk)
{
    using namespace seqan3::detail;

    std::vector<int> in{1, 2, 3, 4, 5, 6};

    chunk_input_iterator it{in};

    it.next_chunk(6);
    auto rng = it.get_chunk();

    EXPECT_EQ(*rng.begin(), in[0]);
    EXPECT_EQ(rng.end(), std::end(in));
}

TEST(direction_iterator_test, test_direction_iterator_chunk_input_iterator_advance_chunk)
{
    using namespace seqan3::detail;

    std::vector<int> in{1, 2, 3, 4, 5, 6};

    chunk_input_iterator it{in};

    it.advance_chunk(3);

    auto rng = it.get_chunk();

    EXPECT_EQ(*rng.begin(), in[3]);
    EXPECT_EQ(rng.end(), std::end(in));
}

TEST(direction_iterator_test, test_direction_iterator_chunk_output_iterator_construction)
{
    using namespace seqan3::detail;

    std::vector<int> out{1, 2, 3, 4, 5, 6};

    EXPECT_NO_THROW(chunk_output_iterator{out});
}

TEST(direction_iterator_test, test_direction_iterator_chunk_output_iterator_increment)
{
    using namespace seqan3::detail;

    std::vector<int> out{1, 2, 3, 4, 5, 6};

    chunk_output_iterator it{out};

    EXPECT_TRUE((std::is_same_v<decltype(++it), chunk_output_iterator<std::vector<int>>&>));
    EXPECT_TRUE((std::is_same_v<decltype(it++), chunk_output_iterator<std::vector<int>>>));
}

TEST(direction_iterator_test, test_direction_iterator_chunk_output_iterator_dereference)
{
    using namespace seqan3::detail;

    std::vector<int> out{1, 2, 3, 4, 5, 6};

    chunk_output_iterator it{out};

    EXPECT_TRUE((std::is_same_v<decltype(*it), chunk_output_iterator<std::vector<int>>&>));
}

TEST(direction_iterator_test, test_direction_iterator_chunk_output_iterator_assign)
{
    using namespace seqan3::detail;

    std::vector<int> out{1, 2, 3, 4, 5, 6};

    chunk_output_iterator it{out};
    *it = 7;

    EXPECT_EQ(out[6], 7);
}

TEST(direction_iterator_test, test_direction_iterator_chunk_output_iterator_next_chunk)
{
    using namespace seqan3::detail;

    std::vector<int> out{1, 2, 3, 4, 5, 6};

    chunk_output_iterator it{out};

    it.next_chunk(6);

    EXPECT_EQ(out.size(), 12u);
    EXPECT_EQ(out[6], 0);
}

TEST(direction_iterator_test, test_direction_iterator_chunk_output_iterator_get_chunk)
{
    using namespace seqan3::detail;

    std::vector<int> out{1, 2, 3, 4, 5, 6};

    chunk_output_iterator it{out};

    auto rng = it.get_chunk();
    EXPECT_EQ(rng.begin(), rng.end());
    EXPECT_EQ(rng.begin(), end(out));

    it.next_chunk(6);
    rng = it.get_chunk();
    EXPECT_NE(rng.begin(), rng.end());
    EXPECT_EQ(rng.begin(), std::begin(out) + 6 );
    EXPECT_EQ(rng.end(), std::end(out));

    std::iota(rng.begin(), rng.end(), 7);

    EXPECT_EQ(out[6], 7);
    EXPECT_EQ(out[7], 8);
    EXPECT_EQ(out[8], 9);
    EXPECT_EQ(out[9], 10);
    EXPECT_EQ(out[10], 11);
    EXPECT_EQ(out[11], 12);

    rng = it.get_chunk();
    EXPECT_NE(rng.begin(), rng.end());
    EXPECT_EQ(rng.begin(), std::begin(out) + 6 );
    EXPECT_EQ(rng.end(), std::end(out));
}

TEST(direction_iterator_test, test_direction_iterator_chunk_output_iterator_advance_chunk)
{
    using namespace seqan3::detail;

    std::vector<int> out{1, 2, 3, 4, 5, 6};

    chunk_output_iterator it{out};

    it.next_chunk(6);

    it.advance_chunk(3);

    auto rng = it.get_chunk();

    EXPECT_NE(rng.begin(), rng.end());
    EXPECT_EQ(rng.begin(), std::begin(out) + 9);
    EXPECT_EQ(rng.end(), std::end(out));
}

TEST(direction_iterator_test, test_direction_iterator_input_iterator_vector)
{
    using namespace seqan3::detail;

    std::vector<char> in{'a', 'c', 'g'};
    auto [r_beg, r_end] = input_iterator(in);
    EXPECT_EQ(*r_beg, 'a');
    EXPECT_NE(r_beg, r_end);
    EXPECT_EQ(r_end, (chunk_input_iterator{in, true}));
}

TEST(direction_iterator_test, test_direction_iterator_output_iterator_vector)
{
    using namespace seqan3::detail;

    std::vector<char> out{'a', 'c', 'g'};

    *output_iterator(out) = 't';
    EXPECT_EQ(out[3], 't');
    EXPECT_EQ(out.size(), 4u);
}
