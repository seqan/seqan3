// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <deque>
#include <list>
#include <map>
#include <vector>

#include <seqan3/core/bit_manipulation.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>

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
    [[maybe_unused]] constexpr aligned_allocator<int, 16> alloc{};
}

TEST(aligned_allocator, conversion_constructor)
{
    aligned_allocator<int, 16> int_alloc{};
    [[maybe_unused]] aligned_allocator<float, 16> float_alloc{int_alloc};
}

TEST(aligned_allocator, request_too_much_memory)
{
    aligned_allocator<int, 16> alloc{};
    EXPECT_THROW((void) alloc.allocate(std::numeric_limits<uint64_t>::max()), std::bad_alloc);
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

    EXPECT_EQ(sizeof(int), 4u);

    EXPECT_EQ(memory_alignment(begin, alignment), 0u);
    EXPECT_EQ(memory_alignment(end,   alignment), 8u);

    EXPECT_EQ(memory_alignment(begin + 1,  alignment), 4u);
    EXPECT_EQ(memory_alignment(begin + 2,  alignment), 8u);
    EXPECT_EQ(memory_alignment(begin + 3,  alignment), 12u);
    EXPECT_EQ(memory_alignment(begin + 4,  alignment), 0u);
    EXPECT_EQ(memory_alignment(begin + 5,  alignment), 4u);
    EXPECT_EQ(memory_alignment(begin + 6,  alignment), 8u);
    EXPECT_EQ(memory_alignment(begin + 7,  alignment), 12u);
    EXPECT_EQ(memory_alignment(begin + 8,  alignment), 0u);
    EXPECT_EQ(memory_alignment(begin + 9,  alignment), 4u);
    EXPECT_EQ(memory_alignment(begin + 10, alignment), 8u);

    alloc.deallocate(begin, size);
}

TEST(aligned_allocator, memory_alignment_bigger_than_default_new_alignment)
{
    size_t size = 10;
    constexpr size_t alignment = detail::next_power_of_two(__STDCPP_DEFAULT_NEW_ALIGNMENT__ + 1);
    aligned_allocator<int, alignment> alloc{};

    int * begin = alloc.allocate(size);
    int * end   = begin + size;

    EXPECT_EQ(sizeof(int), 4u);

    EXPECT_EQ(memory_alignment(begin, alignment), 0u);
    EXPECT_EQ(memory_alignment(end,   alignment), 8u);

    EXPECT_EQ(memory_alignment(begin + 1,  alignment), 4u);
    EXPECT_EQ(memory_alignment(begin + 2,  alignment), 8u);
    EXPECT_EQ(memory_alignment(begin + 3,  alignment), 12u);
    EXPECT_EQ(memory_alignment(begin + 4,  alignment), 16u);
    EXPECT_EQ(memory_alignment(begin + 5,  alignment), 20u);
    EXPECT_EQ(memory_alignment(begin + 6,  alignment), 24u);
    EXPECT_EQ(memory_alignment(begin + 7,  alignment), 28u);
    EXPECT_EQ(memory_alignment(begin + 8,  alignment), 0u);
    EXPECT_EQ(memory_alignment(begin + 9,  alignment), 4u);
    EXPECT_EQ(memory_alignment(begin + 10, alignment), 8u);

    alloc.deallocate(begin, size);
}

struct large_alignment
{
    alignas(64) std::array<int32_t, 2> data{};
};

TEST(aligned_allocator, memory_alignment_with_large_alignment_type)
{
    size_t size = 10;
    aligned_allocator<large_alignment, alignof(large_alignment)> alloc{};

    large_alignment * begin = alloc.allocate(size);
    large_alignment * end   = begin + size;

    constexpr size_t alignment = alignof(large_alignment);
    EXPECT_EQ(sizeof(large_alignment), 64u);
    EXPECT_EQ(alignof(large_alignment), 64u);

    EXPECT_EQ(memory_alignment(begin, alignment), 0u);
    EXPECT_EQ(memory_alignment(end,   alignment), 0u);

    EXPECT_EQ(memory_alignment(begin + 1,  alignment), 0u);
    EXPECT_EQ(memory_alignment(begin + 2,  alignment), 0u);
    EXPECT_EQ(memory_alignment(begin + 3,  alignment), 0u);
    EXPECT_EQ(memory_alignment(begin + 4,  alignment), 0u);
    EXPECT_EQ(memory_alignment(begin + 5,  alignment), 0u);
    EXPECT_EQ(memory_alignment(begin + 6,  alignment), 0u);
    EXPECT_EQ(memory_alignment(begin + 7,  alignment), 0u);
    EXPECT_EQ(memory_alignment(begin + 8,  alignment), 0u);
    EXPECT_EQ(memory_alignment(begin + 9,  alignment), 0u);
    EXPECT_EQ(memory_alignment(begin + 10, alignment), 0u);

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

    EXPECT_EQ(sizeof(int), 4u);

    EXPECT_EQ(memory_alignment(&*begin_it, alignment), 0u);
    EXPECT_EQ(memory_alignment(&*end_it,   alignment), 8u);

    EXPECT_EQ(memory_alignment(&*(++it), alignment), 4u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 8u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 12u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 4u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 8u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 12u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 4u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 8u);
}

TEST(aligned_allocator, in_deque)
{
    size_t size = 10;
    constexpr size_t alignment = 16;
    std::deque<int, aligned_allocator<int, alignment>> container(size);

    auto begin_it = container.begin();
    auto it       = begin_it;
    auto end_it   = container.end();

    EXPECT_EQ(sizeof(int), 4u);

    EXPECT_EQ(memory_alignment(&*begin_it, alignment), 0u);
    EXPECT_EQ(memory_alignment(&*end_it,   alignment), 8u);

    EXPECT_EQ(memory_alignment(&*(++it), alignment), 4u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 8u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 12u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 4u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 8u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 12u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 4u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 8u);
}

TEST(aligned_allocator, in_list)
{
    size_t size = 10;
    constexpr size_t alignment = 16;
    std::list<int, aligned_allocator<int, alignment>> container(size);

    auto begin_it = container.begin();
    auto it       = begin_it;
    auto end_it   = container.end();

    EXPECT_EQ(sizeof(int), 4u);

    EXPECT_EQ(memory_alignment(&*begin_it, alignment), 0u);
    EXPECT_EQ(memory_alignment(&*end_it,   alignment), 0u);

    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
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

    EXPECT_EQ(sizeof(int), 4u);

    EXPECT_EQ(memory_alignment(&*begin_it, alignment), 0u);

    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
    EXPECT_EQ(memory_alignment(&*(++it), alignment), 0u);
}
