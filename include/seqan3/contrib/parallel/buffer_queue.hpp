// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::buffer_queue.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <atomic>
#include <cmath>
#include <mutex>
#include <new>
#include <shared_mutex>
#include <type_traits>
#include <vector>

#include <seqan3/core/parallel/detail/spin_delay.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/concepts>
#include <seqan3/std/new>
#include <seqan3/std/ranges>
#include <seqan3/std/span>

namespace seqan3::contrib
{

//!\cond
enum class queue_op_status : uint8_t
{
    success = 0,
    empty,
    full,
    closed
};

enum struct buffer_queue_policy : uint8_t
{
    fixed,
    dynamic
};

template <std::Semiregular value_t,
          SequenceContainer buffer_t = std::vector<value_t>,
          buffer_queue_policy buffer_policy = buffer_queue_policy::dynamic>
class buffer_queue
{
public:

    using buffer_type     = buffer_t;
    using value_type      = typename buffer_type::value_type;
    using size_type       = typename buffer_type::size_type;
    using reference       = void;
    using const_reference = void;

    // Default constructor sets capacity to 1 (still empty)
    buffer_queue() : buffer_queue{0u}
    {}
    buffer_queue(buffer_queue const &)             = delete;
    buffer_queue(buffer_queue &&)                  = default;
    buffer_queue & operator=(buffer_queue const &) = delete;
    buffer_queue & operator=(buffer_queue &&)      = default;
    ~buffer_queue()                                = default;

    // you can set the initial capacity here
    explicit buffer_queue(size_type const init_capacity)
    {
        data.resize(init_capacity + 1);
        roundSize = static_cast<size_type>(1) << static_cast<size_type>(std::log2(data.size() - 1) + 1);
    }

    template <std::ranges::InputRange range_type>
        requires std::ConvertibleTo<value_type_t<range_type>, value_type>
    buffer_queue(size_type const init_capacity, range_type && r) : buffer_queue{init_capacity}
    {
        std::ranges::copy(r, std::ranges::begin(data));
    }

    /*!\name Waiting operations
     * \{
     */
    template <typename value2_t>
        requires std::ConvertibleTo<value2_t, value_t>
    void push(value2_t && value)
    {
        detail::spin_delay delay{};

        for (;;)
        {
            auto status = try_push(std::forward<value2_t>(value));
            if (status == queue_op_status::closed)
                throw queue_op_status::closed;
            else if (status == queue_op_status::success)
                return;

            assert(status != queue_op_status::empty);
            assert(status == queue_op_status::full);
            delay.wait(); // pause and then try again.
        }
    } // throws if closed

    template <typename value2_t>
        requires std::ConvertibleTo<value2_t, value_t>
    queue_op_status wait_push(value2_t && value)
    {
        detail::spin_delay delay{};

        for (;;)
        {
            auto status = try_push(std::forward<value2_t>(value));
            // wait until queue is not full anymore..
            if (status != queue_op_status::full)
                return status;

            assert(status != queue_op_status::empty);
            assert(status == queue_op_status::full);
            delay.wait(); // pause and then try again.
        }
    }

    value_type value_pop() // throws if closed
    {
        detail::spin_delay delay{};

        value_type value{};
        for (;;)
        {
            auto status = try_pop(value);

            if (status == queue_op_status::closed)
                throw queue_op_status::closed;
            else if (status == queue_op_status::success)
                return value;

            assert(status != queue_op_status::full);
            assert(status == queue_op_status::empty);
            delay.wait(); // pause and then try again.
        }
    }

    queue_op_status wait_pop(value_type & value)
    {
        detail::spin_delay delay{};

        queue_op_status status;
        for (;;)
        {
            status = try_pop(value);

            if (status == queue_op_status::closed || status == queue_op_status::success)
                break;

            assert(status != queue_op_status::full);
            assert(status == queue_op_status::empty);
            delay.wait(); // pause and then try again.
        }
        return status;
    }
    //!\}

    /*!\name Non-waiting operations
     * \{
     */
    template <typename value2_t>
        requires std::ConvertibleTo<value2_t, value_t>
    queue_op_status try_push(value2_t &&);

    queue_op_status try_pop(value_t &);
    //!\}

    /*!\name State operations
     * \{
     */
    void close() noexcept
    {
        closed_flag.store(true, std::memory_order_release);
    }

    bool is_closed() const noexcept
    {
        return closed_flag.load(std::memory_order_acquire);
    }

    bool is_empty() const noexcept
    {
        std::unique_lock write_lock(mutex);
        return headPos == tailPos;
    }

    bool is_full() const noexcept
    {
        std::unique_lock write_lock(mutex);
        return size_impl() == data.max_size();
    }

    size_type size() const noexcept
    {
        std::unique_lock write_lock(mutex);
        return size_impl();
    }
    //!\}
private:

    // Needed in two functions that both acquire exclusive rights.
    size_type size_impl() const noexcept
    {
        size_type mask = roundSize - 1;
        if ((headPos & mask) <= (tailPos & mask))
            return tailPos - headPos;
        else
            return tailPos - headPos - (roundSize - data.size());
    }

    size_type cyclic_increment(size_type value,
                               size_type modulo,
                               size_type roundSize)
    {
        // invariants:
        //   - roundSize is a power of 2
        //   - (value % roundSize) is in [0, modulo)
        //
        // return the next greater value that fulfils the invariants
        // increment write position & roundSize 2^b- 1 -> all bits set.
        if ((++value & (roundSize - 1)) >= modulo)
            value += roundSize - modulo;
        return value;
    }

    template <typename value2_t>
        requires (std::ConvertibleTo<value2_t, value_t>) &&
                 (buffer_policy == buffer_queue_policy::fixed)
    bool overflow(value2_t &&)
    {
        return false;
    }

    template <typename value2_t>
        requires (std::ConvertibleTo<value2_t, value_t>) &&
                 (buffer_policy == buffer_queue_policy::dynamic)
    bool overflow(value2_t && value);

    //!\brief The ring buffer.
    buffer_t data;
    alignas(std::hardware_destructive_interference_size) std::shared_mutex mutable mutex{};
    alignas(std::hardware_destructive_interference_size) std::atomic<size_type>    headPos{0};
    alignas(std::hardware_destructive_interference_size) std::atomic<size_type>    headReadPos{0};
    alignas(std::hardware_destructive_interference_size) std::atomic<size_type>    tailPos{0};
    alignas(std::hardware_destructive_interference_size) std::atomic<size_type>    tailWritePos{0};
    alignas(std::hardware_destructive_interference_size) std::atomic<size_type>    roundSize{0};
    alignas(std::hardware_destructive_interference_size) std::atomic<bool>         closed_flag{false};
};

// Specifies a fixed size buffer queue.
template <std::Semiregular value_t, SequenceContainer buffer_t = std::vector<value_t>>
using fixed_buffer_queue = buffer_queue<value_t, buffer_t, buffer_queue_policy::fixed>;

// Specifies a dynamic size buffer queue (growable).
template <std::Semiregular value_t, SequenceContainer buffer_t = std::vector<value_t>>
using dynamic_buffer_queue = buffer_queue<value_t, buffer_t, buffer_queue_policy::dynamic>;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename value_t, typename buffer_t, buffer_queue_policy buffer_policy>
template <typename value2_t>
    requires (std::ConvertibleTo<value2_t, value_t>) &&
             (buffer_policy == buffer_queue_policy::dynamic)
inline bool buffer_queue<value_t, buffer_t, buffer_policy>::overflow(value2_t && value)
{
    // try to extend capacity
    std::unique_lock write_lock{mutex};
    size_type cap = data.size();
    size_type roundSize = this->roundSize;
    size_type headPos = this->headPos;
    size_type tailPos = this->tailPos;

    assert(tailPos == this->tailWritePos);
    assert(headPos == this->headReadPos);

    bool valueWasAppended = false;

    // did we reach the capacity limit (another thread could have done the upgrade already)?
    if (cyclic_increment(tailPos, cap, roundSize) >= headPos + roundSize)
    {
        if (cap != 0)
        {
            // tailPos & roundSize - 1 == tailPos % capacity
            auto it = std::ranges::begin(data) + (tailPos & (roundSize - 1));
            *it = std::forward<value2_t>(value);
            tailPos = headPos + roundSize;
            valueWasAppended = true;
        }

        assert(tailPos == headPos + roundSize);

        // get positions of head/tail in current data sequence
        size_type headIdx = headPos & (roundSize - 1);
        size_type tailIdx = tailPos & (roundSize - 1);

        // increase capacity
        data.resize(cap + 1);
        size_type delta = data.size() - cap;
        assert(delta == 1);
        roundSize = static_cast<size_type>(1) << ((data.size() > 1) ? static_cast<size_type>(std::log2(data.size() - 1) + 1) : 1);

        std::ranges::move_backward(std::span{data.data() + headIdx, data.data() + cap}, data.data() + data.size());
        if (cap != 0)
        {
            this->headReadPos = this->headPos = headIdx + delta;
            this->tailWritePos = this->tailPos = tailIdx + roundSize;
        }
        this->roundSize = roundSize;
    }
    return valueWasAppended;
}

// ----------------------------------------------------------------------------
// Function try_pop()
// ----------------------------------------------------------------------------

/*
 * @fn ConcurrentQueue#tryPopFront
 * @headerfile <seqan/parallel.h>
 * @brief Try to dequeue a value from a queue.
 *
 * @signature bool tryPopFront(result, queue[, parallelTag]);
 *
 *
 * @param[in,out] queue       A queue.
 * @param[out]    result      The dequeued value (if available).
 * @param[in]     parallelTag The concurrency scheme. If multiple threads dequeue values concurrently this tag must be
 *                            @link ParallelismTags#Parallel @endlink. The more efficient @link ParallelismTags#Serial
 *                            @endlink tag can only be used if one thread calls <tt>popFront</tt> at a time.
 *                            Default is @link ParallelismTags#Parallel @endlink.
 * @return        bool        Returns <tt>true</tt> if a value could be dequeued and <tt>false</tt> otherwise.
 */

//
//  [  ?  ]  [  4  ]  [  3  ]  [  8  ]  [  0  ]  [  x  ]  [  ?  ]
//                       |                          ^
//                       v                          |
//             head            headRead   tail  tailWrite
//
// empty = (head == tail)
// full = (tail + 1 == head)
//
// valid data between  [headRead, tail)
// currently filled    [tail, tailWrite)
// currently removed   [head, headRead)

template <typename value_t, typename buffer_t, buffer_queue_policy buffer_policy>
inline queue_op_status buffer_queue<value_t, buffer_t, buffer_policy>::try_pop(value_t & result)
{
    // try to extract a value
    std::shared_lock  read_lock{mutex};

    size_type cap = data.size();
    size_type roundSize = this->roundSize;
    size_type headReadPos;
    size_type newHeadReadPos;
    detail::spin_delay spinDelay;

    // wait for queue to become filled
    while (true)
    {
        headReadPos = this->headReadPos;
        size_type tailPos = this->tailPos;

        assert(headReadPos <= tailPos);

        // return if queue is empty
        if (headReadPos == tailPos)
        {
            if (is_closed())  // if empty and closed, no more data is expected.
                return queue_op_status::closed;
            return queue_op_status::empty;
        }

        newHeadReadPos = cyclic_increment(headReadPos, cap, roundSize);

        if (this->headReadPos.compare_exchange_weak(headReadPos, newHeadReadPos))
            break;

        spinDelay.wait();
    }

    // extract value and destruct it in the data string
    result = std::ranges::iter_move(data.begin() + (headReadPos & (roundSize - 1)));

    // wait for pending previous reads and synchronize headPos to headReadPos
    {
        detail::spin_delay delay{};
        size_type old = headReadPos;
        while (!this->headPos.compare_exchange_weak(old, newHeadReadPos))
        {
            old = headReadPos;
            delay.wait();  // add adapting delay in case of high contention.
        }
    }

    return queue_op_status::success;
}

// ----------------------------------------------------------------------------
// Function try_push()
// ----------------------------------------------------------------------------

/*
 * @fn ConcurrentQueue#appendValue
 * @headerfile <seqan/parallel.h>
 * @brief Enqueue a value to a queue.
 *
 * @signature void appendValue(queue, val[, expandTag[, parallelTag]);
 *
 *
 * @param[in,out] queue       A queue.
 * @param[in]     val         The value to enqueue.
 * @param[in]     expandTag   The overflow strategy. If @link OverflowStrategyTags#Generous @endlink the queue will be
 *                            automatically resized if the capacity is exceeded, otherwise the thread spinlocks until
 *                            the element can be enqueued.
 *                            Default is the @link DefaultOverflowImplicit @endlink result for the <tt>queue</tt> type.
 * @param[in]     parallelTag The concurrency scheme. If multiple threads enqueue values concurrently this tag must be
 *                            @link ParallelismTags#Parallel @endlink. The more efficient @link ParallelismTags#Serial
 *                            @endlink tag can only be used if one thread calls <tt>appendValue</tt> at a time.
 *                            Default is @link ParallelismTags#Parallel @endlink.
 */
//
template <typename value_t, typename buffer_t, buffer_queue_policy buffer_policy>
template <typename value2_t>
    requires std::ConvertibleTo<value2_t, value_t>
inline queue_op_status buffer_queue<value_t, buffer_t, buffer_policy>::try_push(value2_t && value)
{
    // try to push the value
    {
        detail::spin_delay delay{};

        std::shared_lock read_lock(mutex);

        if (is_closed())
            return queue_op_status::closed;

        size_type cap = data.size();
        size_type roundSize = this->roundSize;

        while (true)
        {
            size_type tailWritePos = this->tailWritePos;
            size_type newTailWritePos = cyclic_increment(tailWritePos, cap, roundSize);
            size_type headPos = this->headPos;

            assert(newTailWritePos <= (headPos + roundSize + 1));

            // break if we have a wrap around, i.e. queue is full
            if (newTailWritePos >= headPos + roundSize)
                break;

            if (this->tailWritePos.compare_exchange_weak(tailWritePos, newTailWritePos))
            {
                auto it = std::ranges::begin(data) + (tailWritePos & (roundSize - 1));
                *it = std::forward<value2_t>(value);
                // here we construct the value into the reserved storage.

                // wait for pending previous writes and synchronise tailPos to tailWritePos
                {
                    detail::spin_delay delay{};
                    size_type old = tailWritePos;
                    while (!this->tailPos.compare_exchange_weak(old, newTailWritePos))
                    {
                        old = tailWritePos;
                        delay.wait();  // add adapting delay in case of high contention.
                    }
                }
                return queue_op_status::success;
            }

            delay.wait();
        }
    }

    // if possible extend capacity and return.
    if (overflow(std::forward<value2_t>(value)))
    {
        return queue_op_status::success;  // always return success, since the queue dynamically resizes and cannot be full.
    }

    // We could not extend the queue so it must be full.
    return queue_op_status::full;
}
//!\endcond
} // namespace seqan3::contrib
