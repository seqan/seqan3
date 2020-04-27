// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <atomic>
#include <chrono>
#include <numeric>
#include <random>
#include <thread>

#include <seqan3/contrib/parallel/buffer_queue.hpp>

template <typename sequential_push_t, typename sequential_pop_t>
void test_buffer_queue_wait_status()
{
    size_t thread_count = std::thread::hardware_concurrency();

    // limit thread count as virtualbox (used by Travis) seems to have problems with thread congestion
    if (thread_count > 4)
        thread_count = 4;

    size_t writer_count = thread_count / 2;
    if constexpr (sequential_push_t::value)
        writer_count = 1;

    if constexpr (sequential_pop_t::value)
        thread_count = writer_count + 1;

    // std::cout << "threads: " << thread_count << ‘\n‘;
    // std::cout << "writers: " << writer_count << ‘\n‘;

    constexpr size_t size_v = 10000;
    seqan3::contrib::dynamic_buffer_queue<uint32_t> queue{100};

    std::atomic<uint32_t> cnt{1};

    // Define the producer section
    auto produce = [&]()
    {
        while (true)
        {
            // Load atomic
            uint32_t intermediate = cnt.fetch_add(1);

            if (intermediate > size_v)
                return;

            // wait semantics
            seqan3::contrib::queue_op_status status = queue.wait_push(intermediate);
            if (status != seqan3::contrib::queue_op_status::success)
                return;
        }
    };

    // Define the consumer section
    std::atomic<uint32_t> sum{0};
    auto consume = [&]() mutable
    {
        uint32_t i = 0;
        while (queue.wait_pop(i) != seqan3::contrib::queue_op_status::closed)
            sum.fetch_add(i, std::memory_order_relaxed);
    };

    // Create producer pool
    std::vector<std::thread> producer_pool;
    for (size_t n = 0; n < writer_count; ++n)
        producer_pool.emplace_back(produce);

    // Create consumer pool
    std::vector<std::thread> consumer_pool;
    for (size_t n = 0; n < thread_count - writer_count; ++n)
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

    EXPECT_EQ(sum.load(), (size_v * (size_v + 1)) / 2);
}

TEST(buffer_queue, spsc_sum)
{
    test_buffer_queue_wait_status<std::true_type, std::true_type>();
}

TEST(buffer_queue, spmc_sum)
{
    test_buffer_queue_wait_status<std::true_type, std::false_type>();
}

TEST(buffer_queue, mpsc_sum)
{
    test_buffer_queue_wait_status<std::false_type, std::true_type>();
}

TEST(buffer_queue, mpmc_sum)
{
    test_buffer_queue_wait_status<std::false_type, std::false_type>();
}

template <typename sequential_push_t, typename sequential_pop_t, seqan3::contrib::buffer_queue_policy buffer_policy>
void test_buffer_queue_wait_throw(size_t initialCapacity)
{
    using queue_t = seqan3::contrib::buffer_queue<size_t, std::vector<size_t>, buffer_policy>;

    queue_t queue{initialCapacity};
    std::vector<size_t> random;
    std::mt19937 rng(0);

    size_t chk_sum = 0;

    random.resize(100000);
    for (size_t i = 0; i < random.size(); ++i)
    {
        random[i] = rng();
        chk_sum ^= random[i];
    }

    volatile std::atomic<size_t> chk_sum2 = 0;
    size_t thread_count = std::thread::hardware_concurrency();

    // limit thread count as virtualbox (used by Travis) seems to have problems with thread congestion
    if (thread_count > 4)
        thread_count = 4;

    size_t writer_count = thread_count / 2;
    if constexpr (sequential_push_t::value)
        writer_count = 1;

    if constexpr (sequential_pop_t::value)
        thread_count = writer_count + 1;

    // std::cout << "threads: " << thread_count << ‘\n‘;
    // std::cout << "writers: " << writer_count << ‘\n‘;

    ASSERT_GE(thread_count, 2u);

    // std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();

    std::vector<std::thread> workers;
    std::atomic<size_t> registered_writer = 0;
    std::atomic<size_t> registered_reader = 0;
    seqan3::contrib::queue_op_status push_status = seqan3::contrib::queue_op_status::success;
    seqan3::contrib::queue_op_status pop_status = seqan3::contrib::queue_op_status::success;
    for (size_t tid = 0; tid < thread_count; ++tid)
    {
        workers.push_back(std::thread([&, tid]()
        {
            // Become writer!
            if (tid < writer_count)
            {
                {  // Wait until all reader are present.
                    seqan3::detail::spin_delay delay{};
                    ++registered_writer;
                    while (registered_reader.load() < (thread_count - writer_count))
                        delay.wait();
                }

                // printf("start writer #%ld\n", tid);
                size_t offset = tid * (random.size() / writer_count);
                size_t offset_end = std::min(static_cast<size_t>((tid + 1) * (random.size() / writer_count)),
                                             random.size());
                for (size_t pos = offset; pos != offset_end; ++pos)
                {
                    try
                    {
                        queue.push(random[pos]);
                    }
                    catch (seqan3::contrib::queue_op_status & ex)
                    {
                        push_status = ex;
                    }
                }
                // printf("stop writer #%ld %lu\n", tid, offset_end - offset);
                // Last writer! No more values will come, so we close the queue.
                ++registered_writer;
                if (registered_writer.load() == (2 * writer_count))
                {
                    queue.close();
                    // printf("writer #%ld closed the queue\n", tid);
                }
            }

            // Become reader!
            if (tid >= writer_count)
            {

                {  // Wait until all writers are setup.
                    seqan3::detail::spin_delay delay{};
                    ++registered_reader;
                    while (registered_writer.load() < writer_count)
                        delay.wait();
                }

                // printf("start reader #%lu\n",  (long unsigned)tid);
                size_t chk_sum_local = 0, cnt = 0;
                for (;;)
                {
                    try
                    {
                        size_t val = queue.value_pop();
                        chk_sum_local ^= val;
                        ++cnt;
                        // if ((cnt & 0xff) == 0)
                        //    printf("%ld ", tid);
                    }
                    catch (seqan3::contrib::queue_op_status & ex)
                    {
                        pop_status = ex;
                        break;
                    }
                }

                chk_sum2.fetch_xor(chk_sum_local);
                // printf("stop reader #%lu %lu\n", static_cast<size_t>(tid), cnt);
            }
        }));
    }

    for (auto & t : workers)
    {
        if (t.joinable())
            t.join();
    }

    // std::chrono::steady_clock::time_point stop = std::chrono::steady_clock::now();
    // double time_span = std::chrono::duration_cast<std::chrono::duration<double> >(stop - start).count();
    // std::cout << "throughput: " << static_cast<size_t>(random.size() / time_span) << " values/s\n";

    EXPECT_EQ(chk_sum, chk_sum2);
    EXPECT_TRUE(push_status == seqan3::contrib::queue_op_status::success);
    EXPECT_TRUE(pop_status == seqan3::contrib::queue_op_status::closed);
}

TEST(buffer_queue, spsc_dynamicsize)
{
    test_buffer_queue_wait_throw<std::true_type, std::true_type, seqan3::contrib::buffer_queue_policy::dynamic>(0u);
}

TEST(buffer_queue, spsc_fixedsize)
{
    test_buffer_queue_wait_throw<std::true_type, std::true_type, seqan3::contrib::buffer_queue_policy::fixed>(30u);
}

TEST(buffer_queue, spmc_dynamicsize)
{
    test_buffer_queue_wait_throw<std::true_type, std::false_type, seqan3::contrib::buffer_queue_policy::dynamic>(0u);
}

TEST(buffer_queue, spmc_fixedsize)
{
    test_buffer_queue_wait_throw<std::true_type, std::false_type, seqan3::contrib::buffer_queue_policy::fixed>(30u);
}

TEST(buffer_queue, mpsc_dynamicsize)
{
    test_buffer_queue_wait_throw<std::false_type, std::true_type, seqan3::contrib::buffer_queue_policy::dynamic>(0u);
}

TEST(buffer_queue, mpsc_fixedsize)
{
    test_buffer_queue_wait_throw<std::false_type, std::true_type, seqan3::contrib::buffer_queue_policy::fixed>(30u);
}

TEST(buffer_queue, mpmc_dynamicsize)
{
    test_buffer_queue_wait_throw<std::false_type, std::false_type, seqan3::contrib::buffer_queue_policy::dynamic>(0u);
}

TEST(buffer_queue, mpmc_fixedsize)
{
    test_buffer_queue_wait_throw<std::false_type, std::false_type, seqan3::contrib::buffer_queue_policy::fixed>(30u);
}
