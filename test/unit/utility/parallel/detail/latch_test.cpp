// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <vector>

#include <seqan3/utility/parallel/detail/latch.hpp>

TEST(latch, arrive_wait)
{
    auto threads = std::thread::hardware_concurrency();
    if (threads > 4)
        threads = 4;

    seqan3::detail::latch completion_latch{threads};
    std::atomic<uint32_t> counter{0};

    auto work = [&]()
    {
        for (unsigned i = 0; i < 1'000'000; ++i)
            ++counter;
        completion_latch.arrive();
    };

    std::vector<std::thread> pool;
    for (unsigned i = 0; i < threads; ++i)
        pool.emplace_back(work);

    completion_latch.wait();

    EXPECT_EQ(counter.load(), 1'000'000 * threads);

    // All threads finished so we can join the threads.
    for (auto & t : pool)
        t.join();
}

TEST(latch, arrive_and_wait)
{
    auto threads = std::thread::hardware_concurrency();
    if (threads > 4)
        threads = 4;

    seqan3::detail::latch completion_latch{threads};
    std::atomic<uint32_t> counter{0};

    auto work = [&]()
    {
        for (unsigned i = 0; i < 1'000'000; ++i)
            ++counter;
        completion_latch.arrive_and_wait();

        EXPECT_EQ(counter.load(), 1'000'000 * threads);
    };

    std::vector<std::thread> pool;
    for (unsigned i = 0; i < threads; ++i)
        pool.emplace_back(work);

    completion_latch.wait();
    EXPECT_EQ(counter.load(), 1'000'000 * threads);
    // All threads finished so we can join the threads.
    for (auto & t : pool)
        t.join();
}
