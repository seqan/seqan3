// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::latch.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <atomic>
#include <cassert>

#include <seqan3/core/parallel/detail/spin_delay.hpp>
#include <seqan3/std/new>

namespace seqan3::detail
{

/*!\brief A single-use synchronisation point to coordinate concurrent threads.
 * \ingroup parallel
 *
 * \details
 *
 * A latch is a thread coordination mechanism that allows any number of threads to block until an expected
 * count is summed (exactly) by threads that arrived at the latch. The expected count is set when the latch is
 * constructed. An individual latch is a single-use object; once the count has been reached, the latch cannot be reused.
 * This implementation uses a lock-free mechanism if the atomic operations on the respective platform are lock-free.
 *
 * \note This adapts the proposal for the c++ standard
 * [P0666R2](http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2018/p0666r2.pdf) for latches and is likely to be
 * changed according to the revisions of the proposal.
 */
class latch
{
public:

    /*!\name Constructors, destructor and assignment
     * \brief Not default constructible nor copyable or movable.
     * \{
     */
    latch()                           = delete;  //!< Deleted.
    latch(latch const &)              = delete;  //!< Deleted.
    latch(latch &&)                   = delete;  //!< Deleted.
    latch & operator=(latch const &)  = delete;  //!< Deleted.
    latch & operator=(latch &&)       = delete;  //!< Deleted.

    //!\brief Destructs the latch and waits for all participating threads to arrive.
    ~latch()
    {
        spin_delay delay{};
        while (num_waiting.load(std::memory_order_acquire) > 0)
            delay.wait();
    }

    /*!\brief Constructs the latch with the expected number of threads.
     * \param expected The number of threads participating in this synchronisation point.
     */
    explicit latch(ptrdiff_t const expected) :
        counter{expected}
    {
        assert(expected >= 0);
        num_waiting.store(0, std::memory_order_relaxed);
    }
    //!\}

    /*!\brief Atomically decrements counter by `n`.
     * \param n The value to subtract.
     *
     * \details
     *
     * Arrives at the synchronisation point with `n` count without waiting. If `n` is less than 0 or greater then the
     * current count, then this operation results in undefined behaviour for all other participating threads.
     *
     * ### Exception
     *
     * Guaranteed not to throw.
     *
     * ### Thread safety
     *
     * Thread-safe.
     */
    void arrive(ptrdiff_t n = 1) noexcept
    {
        assert(counter.load(std::memory_order_acquire) >= n);
        assert(counter.load(std::memory_order_acquire) >= 0);

        counter.fetch_sub(n, std::memory_order_acq_rel);
    }

    /*!\brief Atomically decrements counter by `n` and blocks the calling thread.
     * \param n The value to subtract.
     *
     * \details
     *
     * Arrives at the synchronisation point with `n` count and blocks the calling thread until all participating threads
     * have arrived. If `n` is less than 0 or greater then the current count, then this operation results in undefined
     * behaviour for all other participating threads.
     *
     * ### Exception
     *
     * Guaranteed not to throw.
     *
     * ### Thread safety
     *
     * Thread-safe.
     */
    void arrive_and_wait(ptrdiff_t n = 1) noexcept
    {
        ++num_waiting; // ensure that destructor is not finished in-between the arrive and wait call.
        arrive(n);
        wait();
        --num_waiting;
    }

    /*!\brief Checks if all participating threads have reached the synchronisation point.
     *
     * ### Exception
     *
     * Guaranteed not to throw.
     *
     * ### Thread safety
     *
     * Thread-safe.
     */
    bool try_wait() const noexcept
    {
        return counter.load(std::memory_order_acquire) == 0;
    }

    /*!\brief Waits for all participating threads to arrive at the synchronisation point.
     *
     * \details
     *
     * If counter == 0 returns immediately, otherwise blocks the calling thread at the synchronisation point.
     * This implementation uses a seqan3::spin_delay to wait for the value.
     *
     * ### Exception
     *
     * Guaranteed not to throw.
     *
     * ### Thread safety
     *
     * Thread-safe.
     */
    void wait() const
    {
        ++num_waiting; // register waiting thread to synchronise with destructor.
        spin_delay delay{};
        while (counter.load(std::memory_order_acquire) > 0)
            delay.wait();
        --num_waiting;
    }

private:

    //!\brief The number of participating threads.
    alignas(std::hardware_destructive_interference_size) std::atomic<std::ptrdiff_t>         counter;
    //!\brief The number of waiting threads.
    alignas(std::hardware_destructive_interference_size) mutable std::atomic<std::ptrdiff_t> num_waiting;
};

} // namespace seqan3::detail
