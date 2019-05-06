// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <chrono>
#include <random>
#include <string>

#include <seqan3/contrib/parallel/concurrent_queue.hpp>

using namespace seqan3;

TEST(concurrent_queue, empty)
{
    concurrent_queue<int> queue{};

    EXPECT_TRUE(queue.is_empty());
}

TEST(concurrent_queue, full)
{
    {
        concurrent_queue<int> queue{};
        EXPECT_TRUE(queue.try_push(3) == queue_op_status::full);
    }

    {
        concurrent_queue<int> queue{1};
        EXPECT_TRUE(queue.try_push(3) == queue_op_status::success);
        int x = -1;
        queue.try_pop(x);
    }
}

TEST(concurrent_queue, push_pop)
{
    {
        concurrent_queue<int> queue{};
        int x = -1;
        EXPECT_TRUE(queue.try_push(3) == queue_op_status::full);
        EXPECT_TRUE(queue.try_pop(x) == queue_op_status::empty);
        EXPECT_EQ(x, -1);

        EXPECT_TRUE(queue.try_push(3) == queue_op_status::success);
        EXPECT_TRUE(queue.try_push(6) == queue_op_status::success);

        EXPECT_TRUE(queue.try_pop(x) == queue_op_status::success);
        EXPECT_EQ(x, 3);
        EXPECT_TRUE(queue.try_pop(x) == queue_op_status::success);
        EXPECT_EQ(x, 6);
        EXPECT_TRUE(queue.try_pop(x) == queue_op_status::empty);
        queue.close();
        EXPECT_TRUE(queue.try_pop(x) == queue_op_status::closed);
    }

    for (unsigned i = 0; i < 10; ++i)
    {
        concurrent_queue<int> queue(i);
        int x = -1;
        EXPECT_TRUE(queue.try_pop(x) == queue_op_status::empty);
        EXPECT_TRUE(queue.is_empty());

        queue.try_push(3);
        queue.try_push(6);

        EXPECT_FALSE(queue.is_empty());
        EXPECT_TRUE(queue.try_pop(x) == queue_op_status::success);
        EXPECT_EQ(x, 3);
        EXPECT_TRUE(queue.try_pop(x) == queue_op_status::success);
        EXPECT_EQ(x, 6);
        EXPECT_TRUE(queue.try_pop(x) == queue_op_status::empty);
    }
}

TEST(concurrent_queue, close)
{
    concurrent_queue<int> queue{};
    int x = -1;
    EXPECT_TRUE(queue.try_push(3) == queue_op_status::full);
    EXPECT_TRUE(queue.try_pop(x) == queue_op_status::empty);
    EXPECT_EQ(x, -1);

    EXPECT_TRUE(queue.try_push(3) == queue_op_status::success);
    EXPECT_TRUE(queue.try_push(6) == queue_op_status::success);
    EXPECT_FALSE(queue.is_closed());
    queue.close();
    EXPECT_TRUE(queue.is_closed());
    EXPECT_TRUE(queue.try_pop(x) == queue_op_status::success);
    EXPECT_EQ(x, 3);
    EXPECT_TRUE(queue.try_pop(x) == queue_op_status::success);
    EXPECT_EQ(x, 6);
    EXPECT_TRUE(queue.try_pop(x) == queue_op_status::closed);
}

TEST(concurrent_queue, size)
{
    concurrent_queue<int> queue{};
    int x = -1;
    EXPECT_EQ(queue.size(), 0u);
    EXPECT_TRUE(queue.try_push(3) == queue_op_status::full);
    EXPECT_EQ(queue.size(), 0u);

    EXPECT_TRUE(queue.try_push(3) == queue_op_status::success);
    EXPECT_TRUE(queue.try_push(6) == queue_op_status::success);
    EXPECT_EQ(queue.size(), 2u);
    EXPECT_TRUE(queue.try_pop(x) == queue_op_status::success);
    EXPECT_EQ(x, 3);
    EXPECT_EQ(queue.size(), 1u);
    EXPECT_TRUE(queue.try_pop(x) == queue_op_status::success);
    EXPECT_EQ(x, 6);
    EXPECT_EQ(queue.size(), 0u);
}

TEST(concurrent_queue, non_pod)
{
    // in a queue of capacity 3 try all 3 states of being empty
    for (int ofs = 1; ofs < 10; ++ofs)
    {
        concurrent_queue<std::string> queue{10};

        for (int i = 0; i < ofs; ++i)
        {
            std::string x{};
            std::string y{"al_"};
            y[2] = '0' + i;
            EXPECT_TRUE(queue.try_push(y) == queue_op_status::success);
            EXPECT_TRUE(queue.try_pop(x) == queue_op_status::success);
            EXPECT_EQ(x, y);
        }
        EXPECT_TRUE(queue.is_empty());
    }
}
