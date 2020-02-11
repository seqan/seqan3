// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::views::async_input_buffer.
 */

#pragma once

#include <thread>

#include <seqan3/contrib/parallel/buffer_queue.hpp>
#include <seqan3/range/views/detail.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

//-----------------------------------------------------------------------------
// This is the path a value takes when using this views::
//   urange
// → async_input_buffer_view.buffer               [size n]
// → async_input_buffer_iterator.cached_value     [size 1]
// → user
//-----------------------------------------------------------------------------

namespace seqan3::detail
{

/*!\brief The type returned by seqan3::views::async_input_buffer.
 * \tparam urng_t The underlying range type.
 * \implements std::ranges::InputRange
 * \ingroup views
 */
template <std::ranges::range urng_t>
class async_input_buffer_view : public std::ranges::view_interface<async_input_buffer_view<urng_t>>
{
private:
    static_assert(std::ranges::input_range<urng_t>,
        "The range parameter to async_input_buffer_view must be at least an std::ranges::InputRange.");
    static_assert(std::ranges::view<urng_t>,
        "The range parameter to async_input_buffer_view must model std::ranges::View.");
    static_assert(std::movable<std::ranges::range_value_t<urng_t>>,
        "The range parameter to async_input_buffer_view must have a value_type that is std::Movable.");
    static_assert(std::constructible_from<std::ranges::range_value_t<urng_t>,
                                          std::remove_reference_t<std::ranges::range_reference_t<urng_t>> &&>,
        "The range parameter to async_input_buffer_view must have a value_type that is constructible by a moved "
        "value of its reference type.");

    //!\brief The iterator type for the underlying range.
    using urng_iterator_type = std::ranges::iterator_t<urng_t>;

    //!\brief Buffer and thread and shared between copies of this type.
    struct state
    {
        //!\brief The underlying range.
        urng_t urange;

        //!\brief The buffer queue.
        contrib::fixed_buffer_queue<std::ranges::range_value_t<urng_t>> buffer;

        //!\brief Thread that rebuffers in the background.
        std::thread producer;
    };

    //!\brief Shared holder of the state.
    std::shared_ptr<state> state_ptr = nullptr;

    //!\brief The iterator of the seqan3::detail::async_input_buffer_view.
    class async_input_buffer_iterator;

public:
    /*!\name Constructor, destructor, and assignment.
     * \{
     */
    async_input_buffer_view()                                         = default; //!< Defaulted.
    async_input_buffer_view(async_input_buffer_view const &)                = default; //!< Defaulted.
    async_input_buffer_view(async_input_buffer_view &&)                     = default; //!< Defaulted.
    async_input_buffer_view & operator=(async_input_buffer_view const &)    = default; //!< Defaulted.
    async_input_buffer_view & operator=(async_input_buffer_view &&)         = default; //!< Defaulted.
    ~async_input_buffer_view()                                        = default; //!< Defaulted.

    //!\brief Construction from the underlying view.
    async_input_buffer_view(urng_t _urng, size_t const buffer_size)
    {
        auto deleter = [] (state * p)
        {
            if (p != nullptr)
            {
                p->buffer.close();
                p->producer.join();
                delete p;
            }
        };

        state_ptr = std::shared_ptr<state>(new state{std::move(_urng),
                                                     contrib::fixed_buffer_queue<std::ranges::range_value_t<urng_t>>{buffer_size},
                                                     std::thread{}}, // thread is set/started below, needs rest of state
                                           deleter);

        auto runner = [&state = *state_ptr] ()
        {
            for (auto && val : state.urange)
                if (state.buffer.wait_push(std::move(val)) == contrib::queue_op_status::closed)
                    break;

            state.buffer.close();
        };

        state_ptr->producer = std::thread{runner};
    }

    //!\brief Construction from std::ranges::ViewableRange.
    template <typename other_urng_t>
    //!\cond
    requires !std::same_as<remove_cvref_t<other_urng_t>, async_input_buffer_view> && // prevent recursive instantiation
             std::ranges::viewable_range<other_urng_t> &&
             std::constructible_from<urng_t, ranges::ref_view<std::remove_reference_t<other_urng_t>>>
    //!\endcond
    async_input_buffer_view(other_urng_t && _urng, size_t const buffer_size) :
        async_input_buffer_view{std::views::all(_urng), buffer_size}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */
    /*!\brief Returns an iterator to the current begin of the underlying range.
     *
     * \details
     *
     * ### Thread-Safety
     *
     * It is thread-safe to call this function. Subsequent calls to begin will result in different
     * iterators that are each valid individually. It is thread-safe to operate on different iterators
     * from different threads (however it is not thread-safe to operate on a single iterator from different
     * threads).
     */
    async_input_buffer_iterator begin()
    {
        assert(state_ptr != nullptr);
        return {state_ptr->buffer};
    }

    //!\brief Const-qualified async_input_buffer_view::begin() is deleted, because iterating changes the view.
    async_input_buffer_iterator begin() const = delete;

    //!\copydoc async_input_buffer_view::begin() const
    async_input_buffer_iterator cbegin() const = delete;

    //!\brief Returns a sentinel.
    std::ranges::default_sentinel_t end()
    {
        return std::ranges::default_sentinel;
    }

    //!\brief Const-qualified async_input_buffer_view::end() is deleted, because iterating changes the view.
    std::ranges::default_sentinel_t end() const = delete;

    //!\copydoc async_input_buffer_view::end() const
    std::ranges::default_sentinel_t cend() const = delete;
    //!\}
};

//!\brief The iterator of the seqan3::detail::async_input_buffer_view.
template <typename urng_t>
class async_input_buffer_view<urng_t>::async_input_buffer_iterator
{
    //!\brief The sentinel type to compare to.
    using sentinel_type = std::ranges::default_sentinel_t;

    //!\brief The pointer to the associated view.
    contrib::fixed_buffer_queue<std::ranges::range_value_t<urng_t>> * buffer_ptr = nullptr;

    //!\brief The cached value this iterator holds.
    mutable std::ranges::range_value_t<urng_t> cached_value;

    //!\brief Whether this iterator is at end (the buffer is empty and closed).
    bool at_end = false;

public:

    /*!\name Associated types
    * \{
    */
    //!\brief Difference type.
    using difference_type   = std::iter_difference_t<urng_iterator_type>;
    //!\brief Value type.
    using value_type        = std::iter_value_t<urng_iterator_type>;
    //!\brief Pointer type.
    using pointer           = value_type *;
    //!\brief Reference type.
    using reference         = value_type &;
    //!\brief Iterator category.
    using iterator_category = void;
    //!\brief Iterator concept.
    using iterator_concept  = std::input_iterator_tag;
    //!\}

    /*!\name Construction, destruction and assignment
     * \brief Not explicitly `noexcept` because this depends on construction/copy/... of value_type.
     * \{
     */
    async_input_buffer_iterator()                                                    = default; //!< Defaulted.
    //TODO: delete:
    async_input_buffer_iterator(async_input_buffer_iterator const & rhs)             = default; //!< Defaulted.
    async_input_buffer_iterator(async_input_buffer_iterator && rhs)                  = default; //!< Defaulted.
    //TODO: delete:
    async_input_buffer_iterator & operator=(async_input_buffer_iterator const & rhs) = default; //!< Defaulted.
    async_input_buffer_iterator & operator=(async_input_buffer_iterator && rhs)      = default; //!< Defaulted.
    ~async_input_buffer_iterator()                                          noexcept = default; //!< Defaulted.

    //!\brief Constructing from the underlying seqan3::async_input_buffer_view.
    async_input_buffer_iterator(contrib::fixed_buffer_queue<std::ranges::range_value_t<urng_t>> & buffer) noexcept :
        buffer_ptr{&buffer}
    {
        ++(*this); // cache first value
    }
    //!\}

    /*!\name Access operations
     * \{
     */
    //!\brief Return the cached value.
    reference operator*() const noexcept
    {
        return cached_value;
    }

    //!\brief Returns pointer to the pointed-to object.
    pointer operator->() const noexcept
    {
        return std::addressof(cached_value);
    }
    //!\}

    /*!\name Iterator operations
     * \{
     */
    //!\brief Pre-increment.
    async_input_buffer_iterator & operator++() noexcept
    {
        if (at_end) // TODO unlikely
            return *this;

        assert(buffer_ptr != nullptr);

        if (buffer_ptr->wait_pop(cached_value) == contrib::queue_op_status::closed)
            at_end = true;

        return *this;
    }

    //!\brief Post-increment.
    void operator++(int) noexcept
    {
        ++(*this);
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Compares for equality with sentinel.
    friend constexpr bool operator==(async_input_buffer_iterator const & lhs,
                                     std::ranges::default_sentinel_t const &) noexcept
    {
        return lhs.at_end;
    }

    //!\copydoc operator==
    friend constexpr bool operator==(std::ranges::default_sentinel_t const &,
                                     async_input_buffer_iterator const & rhs) noexcept
    {
        return rhs == std::ranges::default_sentinel_t{};
    }

    //!\brief Compares for inequality with sentinel.
    friend constexpr bool operator!=(async_input_buffer_iterator const & lhs,
                                     std::ranges::default_sentinel_t const &) noexcept
    {
        return !(lhs == std::ranges::default_sentinel_t{});
    }

    //!\copydoc operator!=
    friend constexpr bool operator!=(std::ranges::default_sentinel_t const &,
                                     async_input_buffer_iterator const & rhs) noexcept
    {
        return rhs != std::ranges::default_sentinel_t{};
    }
    //!\}
};

/*!\name Deduction guide.
 * \relates seqan3::detail::async_input_buffer_view
 * \{
 */

//!\brief Deduces the async_input_buffer_view from the underlying range if it is a std::ranges::ViewableRange.
template <std::ranges::viewable_range urng_t>
async_input_buffer_view(urng_t &&, size_t const buffer_size) -> async_input_buffer_view<std::ranges::all_view<urng_t>>;
//!\}

// ============================================================================
//  async_input_buffer_fn (adaptor definition
// ============================================================================

//!\brief Definition of the range adaptor object type for seqan3::views::async_input_buffer.
struct async_input_buffer_fn
{
    //!\brief Store the argument and return a range adaptor closure object.
    constexpr auto operator()(size_t const buffer_size) const
    {
        return detail::adaptor_from_functor{*this, buffer_size};
    }

    /*!\brief Directly return an instance of the view, initialised with the given parameters.
     * \param[in] urange      The underlying range.
     * \param[in] buffer_size The frame that should be used for translation.
     * \returns A range of translated sequence(s).
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, size_t const buffer_size) const
    {
        static_assert(std::ranges::input_range<urng_t>,
            "The range parameter to views::async_input_buffer must be at least an std::ranges::InputRange.");
        static_assert(std::ranges::viewable_range<urng_t>,
            "The range parameter to views::async_input_buffer cannot be a temporary of a non-view range.");
        static_assert(std::movable<std::ranges::range_value_t<urng_t>>,
            "The range parameter to views::async_input_buffer must have a value_type that is std::Movable.");
        static_assert(std::constructible_from<std::ranges::range_value_t<urng_t>,
                                              std::remove_reference_t<std::ranges::range_reference_t<urng_t>> &&>,
            "The range parameter to views::async_input_buffer must have a value_type that is constructible by a moved "
            "value of its reference type.");

        if (buffer_size == 0)
            throw std::invalid_argument{"The buffer_size parameter to views::async_input_buffer must be > 0."};

        return detail::async_input_buffer_view{std::forward<urng_t>(urange), buffer_size};
    }
};

}  // seqan3::detail

//-----------------------------------------------------------------------------
// View shortcut for functor.
//-----------------------------------------------------------------------------

namespace seqan3::views
{
/*!\name General purpose views
 * \{
 */

/*!\brief A view adapter that returns a concurrent-queue-like view over the underlying range.
 * \tparam urng_t         The type of the range being processed. See below for requirements.
 * \param[in,out] urange  The range being processed.
 * \param[in] buffer_size Size of the buffer. Choose the size (> 0) depending on the expected work per element.
 * \returns A view that pre-fetches elements from the underlying range and provides a thread-safe interface.
 *          See below for the properties of the returned range.
 * \ingroup views
 *
 * \details
 *
 * **Header**
 * ```cpp
 * #include <seqan3/range/views/async_input_buffer.hpp>
 * ```
 *
 * ### Summary
 *
 * This view spawns a background thread that pre-fetches elements from the underlying range and stores them in a
 * concurrent queue. Iterating over this view then pops elements out of the queue and returns them.
 * This is primarily useful if dereferencing/incrementing the iterator of the underlying range
 * is expensive, e.g. with SeqAn files which lazily perform I/O.
 *
 * Another advantage of this view is that multiple iterators can be created that are safe to iterate individually,
 * even from different threads, i.e. you can use multiple threads to iterate safely over a single-pass input view
 * with the added benefit of background pre-fetching.
 *
 * In technical terms: this view facilitates a single-producer, multi-consumer design; it's a range interface over
 * a concurrent queue.
 *
 * ### Size of the buffer
 *
 * The `buffer_size` parameter should be chosen depending on the expected work per element, e.g. if the underlying
 * range is an input file over short reads, a buffer size of 100 or 1000 could be beneficial; if on the other hand
 * the file contains genome-sized sequences, it would be better to buffer only a single sequence (buffering 100
 * sequences would result in the entire file being preloaded and likely consuming significant memory).
 *
 * ### Range consumption
 *
 * This view always moves elements from the underlying range into its buffer which means that the elements in
 * the underlying range will be invalidated! For underlying ranges that are single-pass, this makes no difference, but
 * it might be unexpected for multi-pass ranges (std::ranges::forward_range).
 *
 * Typically this adaptor is used when you want to consume the entire underlying range. Destructing
 * this view before all elements have been read will also stop the thread that moves object from the underlying
 * range.
 * **In general, it is not safe to access the underlying range in other contexts once it has been passed
 * to seqan3::views::async_input_buffer.**
 *
 * Note that in addition to the buffer of the view, every iterator has its own one-element-buffer. Dereferencing
 * the iterator returns a reference to the element in the buffer, usually you will want to move this element out
 * of the buffer with std::move std::ranges::iter_move. Incrementing the iterator refills the buffer from the queue
 * inside the view (which in turn is then refilled from the underlying range).
 *
 * ### View properties
 *
 * | concepts and reference type               | `urng_t` (underlying range type)  | `rrng_t` (returned range type)         |
 * |-------------------------------------------|:---------------------------------:|:--------------------------------------:|
 * | std::ranges::input_range                  | *required*                        | *preserved*                            |
 * | std::ranges::forward_range                |                                   | *lost*                                 |
 * | std::ranges::bidirectional_range          |                                   | *lost*                                 |
 * | std::ranges::random_access_range          |                                   | *lost*                                 |
 * | std::ranges::contiguous_range             |                                   | *lost*                                 |
 * |                                           |                                   |                                        |
 * | std::ranges::viewable_range               | *required*                        | *guaranteed*                           |
 * | std::ranges::view                         |                                   | *guaranteed*                           |
 * | std::ranges::sized_range                  |                                   | *lost*                                 |
 * | std::ranges::common_range                 |                                   | *lost*                                 |
 * | std::ranges::output_range                 |                                   | *lost*                                 |
 * | seqan3::const_iterable_range              |                                   | *lost*                                 |
 * |                                           |                                   |                                        |
 * | std::ranges::range_reference_t            |                                   | `std::ranges::range_value_t<urng_t> &` |
 * |                                           |                                   |                                        |
 * | std::iterator_traits \::iterator_category |                                   | *none*                                 |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Thread safety
 *
 * The following operations are **thread-safe**:
 *
 *   * calling `.begin()` and `.end()` on the view returned by this adaptor;
 *   * calling operators on the different iterator objects.
 *
 * Calling operators on the same iterator object from different threads is not safe, i.e. you can pass the view
 * to different threads by reference, and have each of those threads call `begin()` on the view and then
 * perform operations (dereference, increment...) on that iterator from the respective thread; but you
 * cannot call `begin()` in a parent thread, pass the iterator to different threads and operate on that
 * concurrently.
 *
 * ### Example
 *
 * \include test/snippet/range/views/async_input_buffer.cpp
 *
 * Running the snippet could yield the following output:
 *
 * ```
 * Thread: 0x80116bf00     Seq:    seq2
 * Thread: 0x80116bf00     Seq:    seq3
 * Thread: 0x80116ba00     Seq:    seq1
 * Thread: 0x80116bf00     Seq:    seq4
 * Thread: 0x80116bf00     Seq:    seq6
 * Thread: 0x80116ba00     Seq:    seq5
 * Thread: 0x80116bf00     Seq:    seq7
 * Thread: 0x80116ba00     Seq:    seq8
 * Thread: 0x80116bf00     Seq:    seq9
 * Thread: 0x80116bf00     Seq:    seq11
 * Thread: 0x80116bf00     Seq:    seq12
 * Thread: 0x80116ba00     Seq:    seq10
 * ```
 * This shows that indeed elements from the underlying range are processed non-sequentially, that there are two threads
 * and that work is "balanced" between them (one thread processed more element than the other, because its "work"
 * per item happened to be smaller).
 *
 * Note that you might encounter jumbled output if by chance two threads write to the stream at the exact same time.
 *
 * If you remove the line starting with `auto f1 = ...` you will get sequential processing:
 * ```
 * Thread: 0x80116aa00     Seq:    seq1
 * Thread: 0x80116aa00     Seq:    seq2
 * Thread: 0x80116aa00     Seq:    seq3
 * Thread: 0x80116aa00     Seq:    seq4
 * Thread: 0x80116aa00     Seq:    seq5
 * Thread: 0x80116aa00     Seq:    seq6
 * Thread: 0x80116aa00     Seq:    seq7
 * Thread: 0x80116aa00     Seq:    seq8
 * Thread: 0x80116aa00     Seq:    seq9
 * Thread: 0x80116aa00     Seq:    seq10
 * Thread: 0x80116aa00     Seq:    seq11
 * Thread: 0x80116aa00     Seq:    seq12
 * ```
 *
 * Note that even if you have a single processing thread, using this view can still improve performance measurably,
 * because loading of the elements into the buffer (which reads input from disk) happens in a background thread.
 *
 * \hideinitializer
 */
inline constexpr auto async_input_buffer = detail::async_input_buffer_fn{};

//!\}
} // namespace seqan3::views
