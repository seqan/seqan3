// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
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

#include <seqan3/core/bit_manipulation.hpp>
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

// Ringbuffer implementation:
// The underlying buffer has size (number of actual elements + 1). This is a trick to easily check if the queue is empty
// or full. Furthermore, the ring buffer uses 4 pointers. The actual push_back and pop_front position as well as
// the pending push_back and pop_front position. The latter indicates the position that have been advanced concurrently
// by multiple threads from either end.
//
// head: position to read/extract from the queue (first inserted elment) => pop_front_position
// tail: position where to write to/push new elements to the queue => push_back_position
// head_read: The actual position after x threads concurrently popped from the queue. => pending_pop_front_position
// tail_write: The actual position after x threads concurrently pushed to the queue. => pending_push_back_position
//  [  ?  ]  [  4  ]  [  3  ]  [  8  ]  [  0  ]  [  x  ]  [  ?  ]
//                       |                          ^
//                       v                          |
//             head            headRead   tail  tailWrite
//
// valid buffer between  [headRead, tail)
// currently filled      [tail, tailWrite)
// currently removed     [head, headRead)
//
// State: empty = (head == tail)
//  [  ?  ]  [  ?  ]  [  ?  ]  [  ?  ]  [  ?  ]  [  ?  ]  [  ?  ]
//                                        tail
//                                        head
// The head is on the same position as tail.
// This means that currently no element is in the buffer.

// State: full = (tail + 1 == head)
//  [  2  ]  [  4  ]  [  3  ]  [  ?  ]  [  8  ]  [  0  ]  [  7  ]
//                               tail
//                                        head
// The tail is one position before the head.
// This means that currently no element can be added to the buffer since it is full.
// Strategies are to either wait until some elements have been popped or to expand the capacity of the
// queue by one, inserting the element at the current tail position and moving all elements starting from head one
// position to the right.

template <std::semiregular value_t,
          sequence_container buffer_t = std::vector<value_t>,
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
    buffer_queue(buffer_queue &&)                  = delete;
    buffer_queue & operator=(buffer_queue const &) = delete;
    buffer_queue & operator=(buffer_queue &&)      = delete;
    ~buffer_queue()                                = default;

    // you can set the initial capacity here
    explicit buffer_queue(size_type const init_capacity)
    {
        buffer.resize(init_capacity + 1);
        ring_buffer_capacity = seqan3::detail::next_power_of_two(buffer.size());
    }

    template <std::ranges::input_range range_type>
        requires std::convertible_to<std::ranges::range_value_t<range_type>, value_type>
    buffer_queue(size_type const init_capacity, range_type && r) : buffer_queue{init_capacity}
    {
        std::ranges::copy(r, std::ranges::begin(buffer));
    }

    /*!\name Waiting operations
     * \{
     */
    template <typename value2_t>
        requires std::convertible_to<value2_t, value_t>
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
        requires std::convertible_to<value2_t, value_t>
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
        requires std::convertible_to<value2_t, value_t>
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
        return pop_front_position == push_back_position;
    }

    bool is_full() const noexcept
    {
        std::unique_lock write_lock(mutex);
        return is_ring_buffer_exhausted(pop_front_position, push_back_position);
    }

    size_type size() const noexcept
    {
        std::unique_lock write_lock(mutex);
        if (to_buffer_position(pop_front_position) <= to_buffer_position(push_back_position))
        {
            return to_buffer_position(push_back_position) - to_buffer_position(pop_front_position);
        }
        else
        {
            assert(buffer.size() > (to_buffer_position(pop_front_position) - to_buffer_position(push_back_position)));
            return buffer.size() - (to_buffer_position(pop_front_position) - to_buffer_position(push_back_position));
        }
    }
    //!\}
private:

    /*!\brief Checks if the capacity of the ring buffer is exhausted.
     * \param[in] from The thread local position to read from the buffer.
     * \param[in] to The thread local position to write to the buffer.
     *
     * \details
     *
     * The seqan3::contrib::buffer_queue::cyclic_increment ensures that the `push_back_position` is at least
     * `ring_buffer_capacity` many slots ahead of `pop_front_position` if the queue is full.
     */
    constexpr bool is_ring_buffer_exhausted(size_type const from, size_type const to) const
    {
        assert(to <= (from + ring_buffer_capacity + 1)); // The tail cannot overwrite the head.

        return to >= from + ring_buffer_capacity;
    }

    /*!\brief Maps the given position to the respective position within the ring buffer.
     * \param[in] position The position to map to the ring buffer position.
     *
     * \details
     *
     * In an emulated ring buffer the actual position within the buffer can be computed using
     * `position % buffer.size()`. The current implementation uses a trick to compute the modulo by means of a bitwise
     * and-operation which is supposedly faster than modulo. To do so, the position is masked with
     * `ring_buffer_capacity - 1`, where `ring_buffer_capacity` is always the next integer number greater or equal
     * than `buffer.size()` that is also a power-of-two. Accordingly, positions smaller than the size of the buffer
     * won't be changed and positions that reach the buffer size are updated by
     * seqan3::contrib::buffer_queue::cyclic_increment in such a way that the masking still works correctly.
     */
    constexpr size_type to_buffer_position(size_type const position) const
    {
       return position & (ring_buffer_capacity - 1);
    }

    /*!\brief Increments the given position by one and returns the next position in the buffer.
     * \param[in] position The position to increment.
     *
     * \details
     *
     * Increments the given position by one. If the position reached the end of the (linear) buffer it emulates a
     * wrap around by adding `ring_buffer_capacity - buffer.size()` to the current position.
     * With this implementation trick it is possible to compute the modulo to get the respective ring buffer position
     * using the function seqan3::contrib::buffer_queue::to_buffer_position.
     *
     * ### Example
     *
     * Assume the buffer has at the moment space for 12 elements. Accordingly, the ring_buffer_capacity is set to 16.
     * If the next ring buffer write position is 12 we need a wrap around in the ring-buffer in order to start at
     * position 0 of the actual underlying buffer instead. To do so the remaining difference between buffer size and
     * ring buffer capacity is added such that the postion is now a multiple of `ring_buffer_capacity`.
     * Accordingly, the modulo operation to receive the actual ring-buffer position can be replaced by
     * `position & (ring_buffer_capacity - 1)`.
     */
    size_type cyclic_increment(size_type position)
    {
        // invariants:
        //   - ring_buffer_capacity is a power of 2
        //   - (position % ring_buffer_capacity) is in [0, buffer.size())
        //
        // return the next greater position that fulfils the invariants
        if (to_buffer_position(++position) >= buffer.size())
            position += ring_buffer_capacity - buffer.size();  // If the position reached
        return position;
    }

    template <typename value2_t>
        requires (std::convertible_to<value2_t, value_t>) &&
                 (buffer_policy == buffer_queue_policy::fixed)
    bool overflow(value2_t &&)
    {
        return false;
    }

    template <typename value2_t>
        requires (std::convertible_to<value2_t, value_t>) &&
                 (buffer_policy == buffer_queue_policy::dynamic)
    bool overflow(value2_t && value);

    //!\brief The buffer that is used as ring buffer.
    buffer_t buffer;
    alignas(std::hardware_destructive_interference_size) std::shared_mutex mutable mutex{};
    alignas(std::hardware_destructive_interference_size) std::atomic<size_type>    pop_front_position{0};
    alignas(std::hardware_destructive_interference_size) std::atomic<size_type>    pending_pop_front_position{0};
    alignas(std::hardware_destructive_interference_size) std::atomic<size_type>    push_back_position{0};
    alignas(std::hardware_destructive_interference_size) std::atomic<size_type>    pending_push_back_position{0};
    alignas(std::hardware_destructive_interference_size) std::atomic<size_type>    ring_buffer_capacity{0};
    alignas(std::hardware_destructive_interference_size) std::atomic<bool>         closed_flag{false};
};

// Specifies a fixed size buffer queue.
template <std::semiregular value_t, sequence_container buffer_t = std::vector<value_t>>
using fixed_buffer_queue = buffer_queue<value_t, buffer_t, buffer_queue_policy::fixed>;

// Specifies a dynamic size buffer queue (growable).
template <std::semiregular value_t, sequence_container buffer_t = std::vector<value_t>>
using dynamic_buffer_queue = buffer_queue<value_t, buffer_t, buffer_queue_policy::dynamic>;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename value_t, typename buffer_t, buffer_queue_policy buffer_policy>
template <typename value2_t>
    requires (std::convertible_to<value2_t, value_t>) &&
             (buffer_policy == buffer_queue_policy::dynamic)
inline bool buffer_queue<value_t, buffer_t, buffer_policy>::overflow(value2_t && value)
{
    // try to extend capacity
    std::unique_lock write_lock{mutex};

    size_type old_size = buffer.size();
    size_type ring_buffer_capacity = this->ring_buffer_capacity;
    size_type local_front = this->pop_front_position;
    size_type local_back = this->push_back_position;

    // Expects no pending pushes or pops in unique lock.
    assert(local_back == this->pending_push_back_position);
    assert(local_front == this->pending_pop_front_position);

    bool valueWasAppended = false;

    // did we reach the capacity limit (another thread could have done the upgrade already)?
    // buffer is full if tail_pos + 1 == head_pos
    if (is_ring_buffer_exhausted(local_front, cyclic_increment(local_back)))
    {
        // In case of a full queue write the value into the additional slot.
        // Note, that the ring-buffer implementation uses one additional field which is not used except
        // when overflow happens. This invariant is used, to simply check for the full/empty state of the queue.
        if (old_size != 0)
        {
            auto it = std::ranges::begin(buffer) + to_buffer_position(local_back);
            *it = std::forward<value2_t>(value);
            local_back = local_front + ring_buffer_capacity;
            valueWasAppended = true;
        }

        assert(is_ring_buffer_exhausted(local_front, local_back));

        // get positions of head/tail in current buffer sequence
        size_type front_buffer_position = to_buffer_position(local_front);
        size_type back_buffer_position = to_buffer_position(local_back);

        // increase capacity by one and move all elements from current pop_front_position one to the right.
        buffer.resize(old_size + 1);
        ring_buffer_capacity = seqan3::detail::next_power_of_two(buffer.size());
        std::ranges::move_backward(std::span{buffer.data() + front_buffer_position, buffer.data() + old_size},
                                   buffer.data() + buffer.size());

        // Update the pop_front and push_back positions.
        if (old_size != 0)
        {
            this->pending_pop_front_position = this->pop_front_position = front_buffer_position + 1;
            this->pending_push_back_position = this->push_back_position = back_buffer_position + ring_buffer_capacity;
        }
        this->ring_buffer_capacity = ring_buffer_capacity;
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
template <typename value_t, typename buffer_t, buffer_queue_policy buffer_policy>
inline queue_op_status buffer_queue<value_t, buffer_t, buffer_policy>::try_pop(value_t & result)
{
    // try to extract a value
    std::shared_lock  read_lock{mutex};

    size_type local_pending_pop_front_position{};
    size_type next_local_pop_front_position{};
    detail::spin_delay spinDelay{};

    local_pending_pop_front_position = this->pending_pop_front_position;
    // wait for queue to become filled
    while (true)
    {
        size_type local_push_back_position = this->push_back_position;

        assert(local_pending_pop_front_position <= local_push_back_position);

        // Check if queue is empty
        if (local_pending_pop_front_position == local_push_back_position)
        {
            return is_closed() ? queue_op_status::closed : queue_op_status::empty;
        }

        // Get the next ring-buffer position to read from.
        next_local_pop_front_position = cyclic_increment(local_pending_pop_front_position);
        // Did another/other thread(s) already acquired this slot?
        // If yes, try with next position. If not, break and read from aquired position.
        if (this->pending_pop_front_position.compare_exchange_weak(local_pending_pop_front_position,
                                                                   next_local_pop_front_position))
            break;

        spinDelay.wait();
    }

    // Store the value from the aquired read position.
    result = std::ranges::iter_move(buffer.begin() + to_buffer_position(local_pending_pop_front_position));

    // wait for pending previous reads and synchronize pop_front_position to local_pending_pop_front_position
    {
        detail::spin_delay delay{};
        size_type acquired_slot = local_pending_pop_front_position;
        while (!this->pop_front_position.compare_exchange_weak(acquired_slot, next_local_pop_front_position))
        {
            acquired_slot = local_pending_pop_front_position;
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
template <typename value_t, typename buffer_t, buffer_queue_policy buffer_policy>
template <typename value2_t>
    requires std::convertible_to<value2_t, value_t>
inline queue_op_status buffer_queue<value_t, buffer_t, buffer_policy>::try_push(value2_t && value)
{
    // try to push the value
    {
        detail::spin_delay delay{};

        std::shared_lock read_lock(mutex);

        if (is_closed())
            return queue_op_status::closed;

        // Current up to date position to push an element to
        size_type local_pending_push_back_position = this->pending_push_back_position;

        while (true)
        {
            // Get the next potential position to write the value too.
            size_type next_local_push_back_position = cyclic_increment(local_pending_push_back_position);
            size_type local_pop_front_position = this->pop_front_position;

            // Check if there are enough slots to write to.
            // If not either wait or try to overflow if it is a dynamic queue.
            if (is_ring_buffer_exhausted(local_pop_front_position, next_local_push_back_position))
                break;

            // Did another/other thread(s) acquired the current pending position before this thread
            // If yes, try again if not, write into acquired slot.
            if (this->pending_push_back_position.compare_exchange_weak(local_pending_push_back_position,
                                                                       next_local_push_back_position))
            {
                // Current thread acquired the local_pending_push_back_position and can now write the value into the
                // proper slot of the ring buffer.
                auto it = std::ranges::begin(buffer) + to_buffer_position(local_pending_push_back_position);
                *it = std::forward<value2_t>(value);

                // wait for pending previous writes and synchronise push_back_position to
                // local_pending_push_back_position
                {
                    detail::spin_delay delay{};
                    // the slot this thread acquired to write to
                    size_type acquired_slot = local_pending_push_back_position;
                    while (!this->push_back_position.compare_exchange_weak(acquired_slot,
                                                                           next_local_push_back_position))
                    {
                        acquired_slot = local_pending_push_back_position;
                        delay.wait();
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
        return queue_op_status::success; // always return success, since the queue resizes and cannot be full.
    }

    // We could not extend the queue so it must be full.
    return queue_op_status::full;
}
//!\endcond
} // namespace seqan3::contrib
