// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::reader_writer_manager.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <mutex>

#include <seqan3/core/parallel/detail/latch.hpp>
#include <seqan3/core/parallel/detail/spin_delay.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/std/new>

namespace seqan3::detail
{

/*!\brief A strong type to set the writer count of a seqan3::detail::reader_writer_manager.
 * \ingroup parallel
 */
struct writer_count : public detail::strong_type<size_t, writer_count>
{
    //!\brief An alias for the base type.
    using base_t = detail::strong_type<size_t, writer_count>;

    //!\brief Import the base constructor.
    using base_t::base_t;
};

/*!\brief A strong type to set the reader count of a seqan3::detail::reader_writer_manager.
 * \ingroup parallel
 */
struct reader_count : public detail::strong_type<size_t, reader_count>
{
    //!\brief An alias for the base type.
    using base_t = detail::strong_type<size_t, reader_count>;

    //!\brief Import the base constructor.
    using base_t::base_t;
};

/*!\brief A single-use synchronisation point for closable concurrent data structures.
 * \ingroup parallel
 *
 *
 * \details
 *
 * A reader writer manager is a thread coordination mechanism specifically designed to synchronise with a concurrent
 * data structure such as the seqan3::detail::buffer_queue. In particular, the manager is constructed with a
 * seqan3::detail::reader_count and a seqan3::detail::writer_count and offers special functions for readers
 * and writers that arrive at the synchronisation point. If all writers arrived at the synchronisation point
 * the associated concurrent data structure is closed by invoking the member method `close()`. This so called
 * completion phase is triggered  by only one of the participating producer threads. Accordingly, no more data can
 * be added to the concurrent data structure and possibly waiting readers on an empty data structure can be released.
 * The destructor of the reader_writer_manager will wait until all readers and writers arrived at the synchronisation
 * point. Thus, a concurrent data structure can be safely destructed after the reader_writer_manager is destructed.
 *
 * The class offers functions to register participating producer and consumer threads, which will signal if they are
 * done when their destructor is called.
 *
 * The reader_writer_manager class is neither copyable nor movable due to internal members used for inter thread
 * synchronisation. Accordingly, the class is not default constructible.
 */
class reader_writer_manager
{
private:

    //!\brief A strictly scope-based seqan3::detail::reader_writer_manager wrapper for producer threads.
    class [[nodiscard]] scoped_writer_type
    {
    public:
        /*!\name Constructors, destructor and assignment
         * \brief Not default constructible nor copyable or movable.
         * \{
         */
        scoped_writer_type() = delete;                                           //!< Deleted.
        scoped_writer_type(scoped_writer_type const &) = default;                //!< Deleted.
        scoped_writer_type(scoped_writer_type &&) = default;                     //!< Defaulted.
        scoped_writer_type & operator=(scoped_writer_type const &) = default;    //!< Deleted.
        scoped_writer_type & operator=(scoped_writer_type &&) = default;         //!< Defaulted.

        /*!\brief Constructs the scoped writer with the associated manager.
         * \param _manager The seqan3::detail::reader_writer_manager.
         */
        explicit scoped_writer_type(reader_writer_manager & _manager) : manager{_manager}
        {}

        //!\brief Calls writer_arrive on the wrapped latch and destructs.
        ~scoped_writer_type()
        {
            manager.writer_arrive();
        }
        //!\}

        //!\brief The wrapped latch.
        reader_writer_manager & manager;
    };

    //!\brief A strictly scope-based seqan3::detail::reader_writer_manager wrapper for consumer threads.
    class [[nodiscard]] scoped_reader_type
    {
    public:
        /*!\name Constructors, destructor and assignment
         * \{
         */
        scoped_reader_type() = delete;                                           //!< Deleted.
        scoped_reader_type(scoped_reader_type const &) = default;                //!< Deleted.
        scoped_reader_type(scoped_reader_type &&) = default;                     //!< Defaulted.
        scoped_reader_type & operator=(scoped_reader_type const &) = default;    //!< Deleted.
        scoped_reader_type & operator=(scoped_reader_type &&) = default;         //!< Defaulted.

        /*!\brief Constructs the scoped reader with the associated manager.
         * \param _manager The seqan3::detail::reader_writer_manager.
         */
        explicit scoped_reader_type(reader_writer_manager & _manager) : manager{_manager}
        {}

        //!\brief Calls reader_arrive on the wrapped latch and destructs.
        ~scoped_reader_type()
        {
            manager.reader_arrive();
        }
        //!\}

        //!\brief The wrapped latch.
        reader_writer_manager & manager;
    };

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    reader_writer_manager()                                          = delete;  //!< Deleted.
    reader_writer_manager(reader_writer_manager const &)             = delete;  //!< Deleted.
    reader_writer_manager(reader_writer_manager &&)                  = delete;  //!< Deleted.
    reader_writer_manager & operator=(reader_writer_manager const &) = delete;  //!< Deleted.
    reader_writer_manager & operator=(reader_writer_manager &&)      = delete;  //!< Deleted.
    ~reader_writer_manager()                                         = default; //!< Defaulted.

    /*!\brief Constructs the reader_writer_manager with the reader count, writer count and the associated data structure.
     * \tparam concurrent_t The type of the concurrent data structure. Must have a member function `close`.
     * \param rcount The expected number of consumer threads. Must be at least one.
     * \param wcount The expected number of producer threads. Must be at least one.
     * \param ds     The concurrent data structure.
     *
     * \details
     *
     * Initialises a std::contrib::latch for the consumer threads and a std::contrib::latch for the producer threads.
     * Initialises the completion function which calls `ds.close()` on completion. Note, only a reference to the
     * concurrent data structure is stored. The caller must ensure that `ds` is not destructed before this.
     *
     * ### Exception
     *
     * Throws std::invalid_argument if reader count or writer count is less than 1.
     */
    template <typename concurrent_t>
    //!\cond
        requires requires { std::declval<concurrent_t>().close(); }  // requires closable concurrent data structure.
    //!\endcond
    reader_writer_manager(reader_count const rcount, writer_count const wcount, concurrent_t & ds) :
        reader_latch{static_cast<ptrdiff_t>(rcount.get())},
        writer_latch{static_cast<ptrdiff_t>(wcount.get())},
        completion_fn{[&ds] () { ds.close(); }}
    {
        if (rcount.get() < 1 || wcount.get() < 1)
            throw std::invalid_argument{"Both, reader count and writer count must be at least 1."};
    }
    //!\}

    /*!\brief Atomically decrements writer counter by one and blocks the calling thread.
     *
     * \details
     *
     * Arrives at the synchronisation point for producer threads and blocks the calling thread until all participating
     * threads have arrived. If all participating producer threads arrived at the synchronisation point, one of the
     * threads will invoke the completion function.
     *
     * ### Exception
     *
     * Guaranteed not to throw.
     *
     * ### Thread safety
     *
     * Thread-safe.
     */
    void writer_arrive_and_wait() noexcept
    {
        writer_latch.arrive_and_wait();

        std::call_once(flag, completion_fn);
    }

    /*!\brief Atomically decrements writer counter by one.
     *
     * \details
     *
     * Arrives at the synchronisation point for producer threads without blocking the calling thread.
     * If all participating producer threads arrived at the synchronisation point, one of the
     * threads will invoke the completion function.
     *
     * ### Exception
     *
     * Guaranteed not to throw.
     *
     * ### Thread safety
     *
     * Thread-safe.
     */
    void writer_arrive() noexcept
    {
        writer_latch.arrive();

        if (writer_latch.try_wait())
            std::call_once(flag, completion_fn);
    }

    /*!\brief Atomically decrements writer counter by one and blocks the calling thread.
     *
     * \details
     *
     * Arrives at the synchronisation point for consumer threads and blocks the calling thread until all participating
     * threads have arrived.
     *
     * ### Exception
     *
     * Guaranteed not to throw.
     *
     * ### Thread safety
     *
     * Thread-safe.
     */
    void reader_arrive_and_wait() noexcept
    {
        reader_latch.arrive_and_wait();
    }

    /*!\brief Atomically decrements reader counter by one.
     *
     * \details
     *
     * Arrives at the synchronisation point for consumer threads without blocking the calling thread.
     *
     * ### Exception
     *
     * Guaranteed not to throw.
     *
     * ### Thread safety
     *
     * Thread-safe.
     */
    void reader_arrive() noexcept
    {
        reader_latch.arrive();
    }

    /*!\brief Registers the current thread as a producer thread for the monitored resource.
     * \returns A RAII wrapper class that automatically deregisters from the monitored resource when destructed.
     *
     * \details
     *
     * If a thread calls this function it will become a participating producer thread for the monitored concurrent
     * data structure. On destruction of the returned RAII wrapper class
     * seqan3::detail::reader_writer_manager::writer_arrive is called automatically.
     *
     * ### Exception
     *
     * Guaranteed not to throw.
     *
     * ### Thread safety
     *
     * Thread-safe.
     */
    scoped_writer_type register_writer() noexcept
    {
        return scoped_writer_type{*this};
    }

    /*!\brief Registers the current thread as a consumer thread for the monitored resource.
     * \returns A RAII wrapper class that automatically deregisters from the monitored resource when destructed.
     *
     * \details
     *
     * If a thread calls this function it will become a participating consumer thread for the monitored concurrent
     * data structure. On destruction of the returned RAII wrapper class
     * seqan3::detail::reader_writer_manager::reader_arrive is called automatically.
     *
     * ### Exception
     *
     * Guaranteed not to throw.
     *
     * ### Thread safety
     *
     * Thread-safe.
     */
    scoped_reader_type register_reader() noexcept
    {
        return scoped_reader_type{*this};
    }

private:

    //!\brief The internal latch for consumer threads.
    alignas(std::hardware_destructive_interference_size) latch          reader_latch;
    //!\brief The internal latch for producer threads.
    alignas(std::hardware_destructive_interference_size) latch          writer_latch;
    //!\brief This flag ensures that the completion phase is invoked only once.
    alignas(std::hardware_destructive_interference_size) std::once_flag flag;
    //!\brief The stored completion function.
    std::function<void()>                                               completion_fn;
};

} // namespace seqan3::detail
