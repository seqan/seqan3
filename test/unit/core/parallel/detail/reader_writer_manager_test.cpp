// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/contrib/parallel/buffer_queue.hpp>
#include <seqan3/core/parallel/detail/reader_writer_manager.hpp>

TEST(reader_writer_manager, parallel)
{
    size_t threads = std::thread::hardware_concurrency();

    if (threads > 4)
        threads = 4;

    if (threads > 1)
        --threads;  // One less for the parallel reader

    EXPECT_GE(threads, 1u);

    uint64_t job_size = threads * 1000000;

    seqan3::contrib::fixed_buffer_queue<uint32_t> source_queue{job_size};
    seqan3::contrib::fixed_buffer_queue<uint32_t> target_queue{job_size};

    seqan3::detail::reader_writer_manager source_manager{seqan3::detail::reader_count{threads},
                                                         seqan3::detail::writer_count{1},
                                                         source_queue};
    seqan3::detail::reader_writer_manager target_manager{seqan3::detail::reader_count{1},
                                                         seqan3::detail::writer_count{threads},
                                                         target_queue};

    auto work = [&] () mutable
    {
        {
            auto reader_agent = source_manager.register_reader();
            auto writer_agent = target_manager.register_writer();

            for (;;)
            {
                uint32_t val = 0;
                if (source_queue.wait_pop(val) == seqan3::contrib::queue_op_status::closed)
                    return;

                seqan3::contrib::queue_op_status status = target_queue.try_push(val);
                EXPECT_TRUE(status == seqan3::contrib::queue_op_status::success);
            }
        }
    };

    uint64_t counter{0};

    // start the producer of source/consumer of target.
    std::thread t1{[&] ()
    {
        {
            auto writer_agent = source_manager.register_writer();

            // Initialise source_queue.
            for (size_t i = 0; i < job_size; ++i)
            {
                seqan3::contrib::queue_op_status status = source_queue.try_push(i + 1);
                EXPECT_TRUE(status == seqan3::contrib::queue_op_status::success);
            }
            EXPECT_FALSE(source_queue.is_closed());
        }
        EXPECT_TRUE(source_queue.is_closed());

        auto reader_agent = target_manager.register_reader();
        for (;;)
        {
            uint32_t val = 0;
            if (target_queue.wait_pop(val) == seqan3::contrib::queue_op_status::closed)
                return;

            counter += val;
        }
    }};

    // Start the consumer of source/producer of target.
    std::vector<std::thread> pool;
    for (size_t i = 0; i < threads; ++i)
        pool.emplace_back(work);

    t1.join();

    EXPECT_TRUE(target_queue.is_closed());
    EXPECT_TRUE(source_queue.is_closed());
    EXPECT_TRUE(target_queue.is_empty());
    EXPECT_TRUE(source_queue.is_empty());
    EXPECT_EQ(counter, static_cast<uint64_t>((job_size * (job_size + 1)) / 2));

    for (std::thread & t : pool)
        t.join();
}
