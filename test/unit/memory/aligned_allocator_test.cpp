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

#include <seqan3/memory/aligned_allocator.hpp>

#include <deque>
#include <list>
#include <map>
#include <vector>

using namespace seqan3;

// standard construction.
TEST(aligned_allocator, standard_construction)
{
    EXPECT_TRUE((std::is_default_constructible_v<aligned_allocator<int, 16>>));
    EXPECT_TRUE((std::is_trivially_default_constructible_v<aligned_allocator<int, 16>>));
    EXPECT_TRUE((std::is_nothrow_default_constructible_v<aligned_allocator<int, 16>>));
    EXPECT_TRUE((std::is_copy_constructible_v<aligned_allocator<int, 16>>));
    EXPECT_TRUE((std::is_trivially_copy_constructible_v<aligned_allocator<int, 16>>));
    EXPECT_TRUE((std::is_nothrow_copy_constructible_v<aligned_allocator<int, 16>>));
    EXPECT_TRUE((std::is_move_constructible_v<aligned_allocator<int, 16>>));
    EXPECT_TRUE((std::is_trivially_move_constructible_v<aligned_allocator<int, 16>>));
    EXPECT_TRUE((std::is_nothrow_move_constructible_v<aligned_allocator<int, 16>>));
    EXPECT_TRUE((std::is_copy_assignable_v<aligned_allocator<int, 16>>));
    EXPECT_TRUE((std::is_trivially_copy_assignable_v<aligned_allocator<int, 16>>));
    EXPECT_TRUE((std::is_nothrow_copy_assignable_v<aligned_allocator<int, 16>>));
    EXPECT_TRUE((std::is_move_assignable_v<aligned_allocator<int, 16>>));
    EXPECT_TRUE((std::is_trivially_move_assignable_v<aligned_allocator<int, 16>>));
    EXPECT_TRUE((std::is_nothrow_move_assignable_v<aligned_allocator<int, 16>>));
}

TEST(aligned_allocator, constexpr_constructor)
{
    constexpr aligned_allocator<int, 16> alloc{};
}

TEST(aligned_allocator, conversion_constructor)
{
    aligned_allocator<int, 16> int_alloc{};
    aligned_allocator<float, 16> float_alloc{int_alloc};
}

size_t memory_alignment(void * value, size_t alignment)
{
   return (((size_t)value) & (alignment-1));
}

TEST(aligned_allocator, memory_alignment)
{
    size_t size = 10;
    constexpr size_t alignment = 16;
    aligned_allocator<int, alignment> alloc{};

    int * begin = alloc.allocate(size);
    int * end   = begin + size;

    EXPECT_EQ(sizeof(int), 4);

    EXPECT_EQ(memory_alignment(begin, alignment), 0);
    EXPECT_EQ(memory_alignment(end,   alignment), 8);

    EXPECT_EQ(memory_alignment(begin + 1,  alignment), 4);
    EXPECT_EQ(memory_alignment(begin + 2,  alignment), 8);
    EXPECT_EQ(memory_alignment(begin + 3,  alignment), 12);
    EXPECT_EQ(memory_alignment(begin + 4,  alignment), 0);
    EXPECT_EQ(memory_alignment(begin + 5,  alignment), 4);
    EXPECT_EQ(memory_alignment(begin + 6,  alignment), 8);
    EXPECT_EQ(memory_alignment(begin + 7,  alignment), 12);
    EXPECT_EQ(memory_alignment(begin + 8,  alignment), 0);
    EXPECT_EQ(memory_alignment(begin + 9,  alignment), 4);
    EXPECT_EQ(memory_alignment(begin + 10, alignment), 8);

    alloc.deallocate(begin, size);
}

TEST(aligned_allocator, in_vector)
{
    size_t size = 10;
    constexpr size_t alignment = 16;
    std::vector<int, aligned_allocator<int, alignment>> container(size);

    auto begin_it = container.begin();
    auto it       = begin_it;
    auto end_it   = container.end();

    EXPECT_EQ(sizeof(int), 4);

    EXPECT_EQ(memory_alignment(&*begin_it, alignment), 0);
    EXPECT_EQ(memory_alignment(&*end_it,   alignment), 8);

    EXPECT_EQ(memory_alignment(&*(++it), alignment), 4);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 8);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 12);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 4);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 8);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 12);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 4);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 8);
}

TEST(aligned_allocator, in_deque)
{
    size_t size = 10;
    constexpr size_t alignment = 16;
    std::deque<int, aligned_allocator<int, alignment>> container(size);

    auto begin_it = container.begin();
    auto it       = begin_it;
    auto end_it   = container.end();

    EXPECT_EQ(sizeof(int), 4);

    EXPECT_EQ(memory_alignment(&*begin_it, alignment), 0);
    EXPECT_EQ(memory_alignment(&*end_it,   alignment), 8);

    EXPECT_EQ(memory_alignment(&*(++it), alignment), 4);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 8);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 12);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 4);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 8);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 12);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 4);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 8);
}

TEST(aligned_allocator, in_list)
{
    size_t size = 10;
    constexpr size_t alignment = 16;
    std::list<int, aligned_allocator<int, alignment>> container(size);

    auto begin_it = container.begin();
    auto it       = begin_it;
    auto end_it   = container.end();

    EXPECT_EQ(sizeof(int), 4);

    EXPECT_EQ(memory_alignment(&*begin_it, alignment), 0);
    EXPECT_EQ(memory_alignment(&*end_it,   alignment), 0);

    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
}

TEST(aligned_allocator, in_map)
{
    constexpr size_t alignment = 16;
    using key_type = char;
    using value_type = int;
    using allocator = aligned_allocator<std::pair<const key_type, value_type>, alignment>;
    std::map<key_type, value_type, std::less<key_type>, allocator> container
    {
        {0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}, {6, 6}, {7, 7}, {8, 8}, {9, 9}
    };

    auto begin_it = container.begin();
    auto it       = begin_it;

    EXPECT_EQ(sizeof(int), 4);

    EXPECT_EQ(memory_alignment(&*begin_it, alignment), 0);

    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0);
}
