// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::execution_handler_parallel.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <functional>
#include <ranges>
#include <thread>
#include <type_traits>
#include <vector>

#include <seqan3/contrib/parallel/buffer_queue.hpp>
#include <seqan3/utility/parallel/detail/reader_writer_manager.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

namespace seqan3::detail
{

/*!\brief Handles the parallel execution of algorithms.
 * \ingroup core_algorithm
 *
 * \details
 *
 * This execution handler implements a non-blocking execute method. This means a call to
 * seqan3::detail::execution_handler_parallel::execute will invoke the algorithm asynchronously.
 * This handler can be used in combination with the seqan3::detail::algorithm_executor_blocking to invoke the
 * algorithms on the given algorithm input.
 *
 * ### Concurrency
 *
 * This class maintains a thread pool and a concurrent queue to execute the algorithm tasks asynchronously.
 * On construction the active consumer threads are spawned in the thread pool and concurrently start fetching
 * algorithm tasks from the concurrent queue. At the same time only one producer thread is allowed to asynchronously
 * submit new algorithm tasks.
 *
 * \note Instances of this class are not copyable.
 *
 * \warning This class is only thread-safe in a single producer context. Multiple consumers are allowed.
 *          Concurrent invocation of the interfaces are undefined behaviour.
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
            state->thread_pool.emplace_back(
                [q]()
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

    /*!\brief Constructs the execution handler spawning 1 thread.
     *
     * \details
     *
     * ### Why only 1 thread?
     *
     * This class is not public. It handles the thread pool when, e.g., using the alignment or search algorithms in
     * parallel via the config. This config requires a value (no default), hence the number of threads is always
     * set by the user.
     *
     * When we use an algorithm in parallel, we also default construct a execution_handler_parallel along the way. If
     * the default is set to use all threads, we have to generate the thread pool and a queue. However, this default
     * constructed execution_handler_parallel is immediately moved away and destructed.
     */
    execution_handler_parallel() : execution_handler_parallel{1u}
    {}

    execution_handler_parallel(execution_handler_parallel const &) = delete;             //!< Deleted.
    execution_handler_parallel(execution_handler_parallel &&) = default;                 //!< Defaulted.
    execution_handler_parallel & operator=(execution_handler_parallel const &) = delete; //!< Deleted.
    execution_handler_parallel & operator=(execution_handler_parallel &&) = default;     //!< Defaulted.
    ~execution_handler_parallel() = default;                                             //!< Defaulted.

    //!\}

    /*!\brief Asynchronously schedules a new algorithm task with the given input and callback.
     * \tparam algorithm_t The type of the algorithm; must model std::copy_constructible and std::invocable with
     *                     the given input type as first argument and the callback type as second argument.
     * \tparam algorithm_input_t The input type to invoke the algorithm with (see below for requirements on this type).
     * \tparam callback_t The type of the callable invoked by the algorithm after generating a new result; must model
     *                    std::copy_constructible.
     *
     * \param[in] algorithm The algorithm to invoke.
     * \param[in] input The input of the algorithm.
     * \param[in] callback A callable which will be invoked on each result generated by the algorithm.
     *
     * \details
     *
     * Inside the function the algorithm and the callback are captured as copies to the sate of a lambda function
     * which wraps the task that is stored on the concurrent queue and asynchronously executed. The algorithm input
     * type, however, is perfectly forwarded if `input` is a lvalue-reference or moved if it is a rvalue-reference.
     * Accordingly, the `algorithm_input_t` must either be a lvalue_reference or std::move_constructible.
     */
    template <std::copy_constructible algorithm_t, typename algorithm_input_t, std::copy_constructible callback_t>
        requires std::invocable<algorithm_t, algorithm_input_t, callback_t>
              && (std::is_lvalue_reference_v<algorithm_input_t> || std::move_constructible<algorithm_input_t>)
    void execute(algorithm_t && algorithm, algorithm_input_t && input, callback_t && callback)
    {
        assert(state != nullptr);

        // Note: Unfortunately, we can't use std::forward_as_tuple here because a std::function object (`task_type`)
        // cannot be constructed if the tuple element type is a rvalue-reference.
        // So we capture the input as a `tuple<algorithm_input_t>` which either is a lvalue reference or has no
        // reference type according to the reference collapsing rules of forwarding references.
        // Then we forward the input into the tuple which either just stores the reference or the input is moved into
        // the tuple. When the task is executed by some thread the stored input will either be forwarded as a
        // lvalue-reference to the algorithm or the input is moved into the algorithm from the tuple. This is valid
        // since the task is executed only once by the parallel execution handler.
        // Here is a discussion about the problem on stackoverflow:
        // https://stackoverflow.com/questions/26831382/capturing-perfectly-forwarded-variable-in-lambda/

        // Asynchronously pushes the algorithm job as a task to the queue.
        // Note: that lambda is mutable, s.t. we can move out the content of input_tpl
        task_type task =
            [=, input_tpl = std::tuple<algorithm_input_t>{std::forward<algorithm_input_t>(input)}]() mutable
        {
            using forward_input_t = std::tuple_element_t<0, decltype(input_tpl)>;
            algorithm(std::forward<forward_input_t>(std::get<0>(input_tpl)), std::move(callback));
        };

        [[maybe_unused]] contrib::queue_op_status status = state->queue.wait_push(std::move(task));
        assert(status == contrib::queue_op_status::success);
    }

    /*!\brief Asynchronously executes the algorithm for every element of the given input range.
     * \tparam algorithm_t The type of the algorithm.
     * \tparam algorithm_input_range_t The input range type.
     * \tparam callback_t The type of the callable invoked by the algorithm after generating a new result.
     *
     * \param[in] algorithm The algorithm to invoke.
     * \param[in] input_range The input range to process asynchronously.
     * \param[in] callback A callable which will be invoked on each result generated by the algorithm for a given input.
     *
     * \details
     *
     * Effectively calls seqan3::detail::execution_handler_parallel::execute on every element of the given input
     * range. For every element, a work task is generated and queued for processing by the threads spawned at the
     * construction of this execution handler.
     * The call blocks until all elements have been processed.
     */
    template <std::copy_constructible algorithm_t,
              std::ranges::input_range algorithm_input_range_t,
              std::copy_constructible callback_t>
        requires std::invocable<algorithm_t, std::ranges::range_reference_t<algorithm_input_range_t>, callback_t>
    void bulk_execute(algorithm_t && algorithm, algorithm_input_range_t && input_range, callback_t && callback)
    {
        for (auto && input : input_range)
            execute(algorithm, std::forward<decltype(input)>(input), callback);

        wait();
    }

    //!\brief Waits until all submitted algorithm jobs have been completed.
    void wait()
    {
        assert(state != nullptr);

        state->stop_and_wait();
    }

private:
    /*!\brief An internal state stored on the heap to allow safe move construction/assignment of the class.
     *
     * \details
     *
     * ### Thread safety
     *
     * This class is only intended for use with a single producer model.
     */
    class internal_state
    {
    public:
        /*!\name Constructors, destructor and assignment
        * \brief Instances of this class are not copyable and not movable.
        * \{
        */
        internal_state() = default;                                  //!< Defaulted.
        internal_state(internal_state const &) = delete;             //!< Deleted.
        internal_state(internal_state &&) = delete;                  //!< Deleted.
        internal_state & operator=(internal_state const &) = delete; //!< Deleted.
        internal_state & operator=(internal_state &&) = delete;      //!< Deleted.

        //!\brief Waits for threads to finish.
        ~internal_state()
        {
            stop_and_wait();
        }
        //!\}

        /*!\brief Waits until all threads have been joined.
         *
         * \details
         *
         * ### Thread safety
         *
         * This function is not thread-safe.
         */
        void stop_and_wait()
        {
            queue.close();

            for (auto & t : thread_pool)
            {
                if (t.joinable())
                    t.join();
            }
        }

        //!\brief The thread pool.
        std::vector<std::thread> thread_pool{};
        //!\brief The concurrent queue containing the algorithms to process.
        contrib::fixed_buffer_queue<task_type> queue{10000};
    };

    //!\brief Manages the internal state.
    std::unique_ptr<internal_state> state{nullptr};
};

} // namespace seqan3::detail
