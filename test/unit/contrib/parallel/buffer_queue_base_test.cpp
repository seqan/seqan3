// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <chrono>
#include <random>
#include <string>

#include <seqan3/contrib/parallel/buffer_queue.hpp>

TEST(buffer_queue, empty)
{
    seqan3::contrib::dynamic_buffer_queue<int> queue{};

    EXPECT_TRUE(queue.is_empty());
}

TEST(buffer_queue, full)
{
    {
        seqan3::contrib::dynamic_buffer_queue<int> queue{};
        EXPECT_TRUE(queue.try_push(3) == seqan3::contrib::queue_op_status::success);
        int x = -1;
        EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::success);
    }

    {
        seqan3::contrib::fixed_buffer_queue<int> queue{};
        EXPECT_TRUE(queue.try_push(3) == seqan3::contrib::queue_op_status::full);
    }

    {
        seqan3::contrib::fixed_buffer_queue<int> queue{2};
        EXPECT_TRUE(queue.try_push(3) == seqan3::contrib::queue_op_status::success);
        EXPECT_TRUE(queue.try_push(6) == seqan3::contrib::queue_op_status::success);
        EXPECT_TRUE(queue.try_push(9) == seqan3::contrib::queue_op_status::full);
        int x = -1;
        EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::success);
        EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::success);
    }
}

TEST(buffer_queue, push_pop)
{
    {
        seqan3::contrib::dynamic_buffer_queue<int> queue{};
        int x = -1;
        EXPECT_TRUE(queue.try_push(3) == seqan3::contrib::queue_op_status::success);
        EXPECT_TRUE(queue.try_push(6) == seqan3::contrib::queue_op_status::success);

        EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::success);
        EXPECT_EQ(x, 3);
        EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::success);
        EXPECT_EQ(x, 6);
        EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::empty);
        queue.close();
        EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::closed);
    }

    for (unsigned i = 0; i < 10; ++i)
    {
        seqan3::contrib::dynamic_buffer_queue<int> queue(i);
        int x = -1;
        EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::empty);
        EXPECT_TRUE(queue.is_empty());

        queue.try_push(3);
        queue.try_push(6);

        EXPECT_FALSE(queue.is_empty());
        EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::success);
        EXPECT_EQ(x, 3);
        EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::success);
        EXPECT_EQ(x, 6);
        EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::empty);
    }
}

TEST(buffer_queue, close)
{
    seqan3::contrib::dynamic_buffer_queue<int> queue{};
    int x = -1;

    EXPECT_TRUE(queue.try_push(3) == seqan3::contrib::queue_op_status::success);
    EXPECT_TRUE(queue.try_push(6) == seqan3::contrib::queue_op_status::success);
    EXPECT_FALSE(queue.is_closed());
    queue.close();
    EXPECT_TRUE(queue.is_closed());
    EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::success);
    EXPECT_EQ(x, 3);
    EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::success);
    EXPECT_EQ(x, 6);
    EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::closed);
}

TEST(buffer_queue, size)
{
    seqan3::contrib::dynamic_buffer_queue<int> queue{};
    int x = -1;
    EXPECT_EQ(queue.size(), 0u);

    EXPECT_TRUE(queue.try_push(3) == seqan3::contrib::queue_op_status::success);
    EXPECT_TRUE(queue.try_push(6) == seqan3::contrib::queue_op_status::success);
    EXPECT_EQ(queue.size(), 2u);
    EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::success);
    EXPECT_EQ(x, 3);
    EXPECT_EQ(queue.size(), 1u);
    EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::success);
    EXPECT_EQ(x, 6);
    EXPECT_EQ(queue.size(), 0u);
}

TEST(buffer_queue, non_pod)
{
    // in a queue of capacity 3 try all 3 states of being empty
    for (int ofs = 1; ofs < 10; ++ofs)
    {
        seqan3::contrib::fixed_buffer_queue<std::string> queue{10};

        for (int i = 0; i < ofs; ++i)
        {
            std::string x{};
            std::string y{"al_"};
            y[2] = '0' + i;
            EXPECT_TRUE(queue.try_push(y) == seqan3::contrib::queue_op_status::success);
            EXPECT_TRUE(queue.try_pop(x) == seqan3::contrib::queue_op_status::success);
            EXPECT_EQ(x, y);
        }
        EXPECT_TRUE(queue.is_empty());
    }
}
