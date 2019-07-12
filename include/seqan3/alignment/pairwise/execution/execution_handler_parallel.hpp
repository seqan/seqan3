// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::execution_handler_parallel.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <thread>
#include <type_traits>
#include <vector>

#include <seqan3/contrib/parallel/buffer_queue.hpp>
#include <seqan3/core/parallel/detail/reader_writer_manager.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief Handles the parallel execution of alignments.
 * \ingroup execution
 *
 * \details
 *
 * This class maintains a thread pool and a concurrent queue to manage the alignment tasks.
 * On construction the active consumer threads are spawned in the thread pool and concurrently start fetching
 * alignment tasks from the concurrent queue. At the same time only one producer thread is allowed to asynchronously
 * submit new alignment tasks.
 *
 * \note Instances of this class are not copyable.
 *
 * \attention This class cannot be reused for multiple calls. For this to work, it requires barriers and a queue that
 *            can be reopened.
 */
class execution_handler_parallel
{
private:
    //!\brief The type erased task type.
    using task_type = std::function<void()>;

public:
    /*!\name Constructors, destructor and assignment
     * \brief Instances of this class are not copyable.
     * \{
     */

    /*!\brief Constructs the execution handler spawning `thread_count` many threads.
     * \param thread_count The number of threads to spawn.
     *
     * \details
     *
     * Spawns `thread_count` many threads processing the tasks in the queue in parallel.
     */
    execution_handler_parallel(size_t const thread_count) : state{std::make_unique<internal_state>()}
    {
        auto * q = &(state->queue);
        for (size_t i = 0; i < thread_count; ++i)
        {
            state->thread_pool.emplace_back([q] ()
            {
                for (;;)
                {
                    task_type task;
                    if (q->wait_pop(task) == contrib::queue_op_status::closed)
                        return;

                    task();
                }
            });
        }
    }

    //!\brief Constructs the execution handler spawning std::thread::hardware_concurrency many threads.
    execution_handler_parallel() : execution_handler_parallel{std::thread::hardware_concurrency()}
    {}

    execution_handler_parallel(execution_handler_parallel const &) = delete;                 //!< Deleted.
    execution_handler_parallel(execution_handler_parallel &&) = default;                     //!< Defaulted.
    execution_handler_parallel & operator=(execution_handler_parallel const &) = delete;     //!< Deleted.
    execution_handler_parallel & operator=(execution_handler_parallel &&) = default;         //!< Defaulted.

    //!\brief Waits for threads to finish.
    ~execution_handler_parallel()
    {
        if (state != nullptr)
            wait();
    }
    //!\}

    /*!\brief Invokes the passed alignment instance in a non-blocking manner.
     * \copydetails execution_handler_sequential::execute
     */
    template <typename fn_type, typename first_range_type, typename second_range_type, typename delegate_type>
    //!\cond
        requires std::Invocable<fn_type, size_t const, first_range_type, second_range_type> &&
                 std::Invocable<delegate_type, std::invoke_result_t<fn_type,
                                                                    size_t const,
                                                                    first_range_type,
                                                                    second_range_type>>
    //!\endcond
    void execute(fn_type && func,
                 size_t const idx,
                 first_range_type first_range,
                 second_range_type second_range,
                 delegate_type && delegate)
    {
        static_assert(std::ranges::View<first_range_type>, "Expected a view!");
        static_assert(std::ranges::View<second_range_type>, "Expected a view!");

        assert(state != nullptr);

        task_type task = [=, first_range = std::move(first_range), second_range = std::move(second_range)] ()
        {
            delegate(func(idx, std::move(first_range), std::move(second_range)));
        };
        // Asynchronously pushes the task to the queue.
        [[maybe_unused]] contrib::queue_op_status status = state->queue.wait_push(std::move(task));
        assert(status == contrib::queue_op_status::success);
    }

    //!\brief Waits until all submitted alignment jobs have been processed.
    void wait()
    {
        assert(state != nullptr);

        if (!state->is_waiting)
        {
            state->is_waiting = true;
            state->queue.close();

            for (auto & t : state->thread_pool)
            {
                if (t.joinable())
                    t.join();
            }
        }
    }

private:
    //!\brief An internal state stored on the heap to allow safe move construction/assignment of the class.
    struct internal_state
    {
        //!\brief The thread pool.
        std::vector<std::thread>                 thread_pool{};
        //!\brief The concurrent queue containing the alignments to process.
        contrib::fixed_buffer_queue<task_type>   queue{10000};
        //!\brief Flag to check if wait was already invoked.
        bool                                     is_waiting{false};
    };

    //!\brief Manages the internal state.
    std::unique_ptr<internal_state> state{nullptr};
};

} // namespace seqan3
