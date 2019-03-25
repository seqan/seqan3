// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <atomic>
#include <numeric>
#include <thread>

#include <seqan3/contrib/conqueue/buffer_queue.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

void test_buffer_queue(size_t num_producer, size_t num_consumer)
{
    constexpr size_t size_v = 10000;
    contrib::buffer_queue<uint32_t> queue{100};

    std::atomic<uint32_t> cnt{0};

    // Define the producer section
    auto produce = [&]()
    {
        while (true)
        {
            // Load atomic
            uint32_t intermediate = cnt.fetch_add(1);

            if (intermediate >= size_v)
                return;

            // wait semantics
            contrib::queue_op_status status = queue.wait_push(intermediate);
            if (status != contrib::queue_op_status::success)
                return;
        }
    };

    // Define the consumer section
    std::atomic<uint32_t> sum{0};
    auto consume = [&]() mutable
    {
        uint32_t i = 0;
        while (queue.wait_pop(i) != contrib::queue_op_status::closed)
            sum.fetch_add(i, std::memory_order_relaxed);
    };

    // Create producer pool
    std::vector<std::thread> producer_pool;
    for (size_t n = 0; n < num_producer; ++n)
        producer_pool.emplace_back(produce);

    // Create consumer pool
    std::vector<std::thread> consumer_pool;
    for (size_t n = 0; n < num_consumer; ++n)
        consumer_pool.emplace_back(consume);


    for (auto & t : producer_pool)
    {
        if (t.joinable())
            t.join();
    }
    // Notify queue that no more work is going to be added.
    queue.close();

    for (auto & t : consumer_pool)
    {
        if (t.joinable())
            t.join();
    }

    auto v = std::view::iota(0, size_v) | std::view::common;
    EXPECT_EQ(sum.load(), (std::accumulate(v.begin(), v.end(), 0u)));
}

TEST(buffer_queue, single_producer_single_consumer)
{
    test_buffer_queue(1, 1);
}

TEST(buffer_queue, single_producer_multiple_consumer)
{
    test_buffer_queue(1, 4);
}

TEST(buffer_queue, multiple_producer_single_consumer)
{
    test_buffer_queue(4, 1);
}

TEST(buffer_queue, multiple_producer_multiple_consumer)
{
    test_buffer_queue(4, 4);
}
