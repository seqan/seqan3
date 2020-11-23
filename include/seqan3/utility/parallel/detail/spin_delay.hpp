// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::spin_delay.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <thread>

#if defined(__SSE2__)  // AMD and Intel
#include <xmmintrin.h>  // _mm_pause()
#endif

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{
/*!\brief A delay for threads waiting for a shared resource.
 * \ingroup parallel
 *
 * \details
 *
 * This spin delay is used in combination with spin locks. If the thread is waiting for a shared resource
 * it will usually waste CPU cycles until it can acquire the lock. Especially if there is a high contention on the
 * lock this can become a performance bottleneck. This spin delay allows for more efficient
 * spinning phases using a hybrid spinning approach. At first it does active spinning until a certain number of
 * wait cycles has been reached. After this the thread yields.
 */
class spin_delay
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr spin_delay()                               noexcept = default;  //!< Defaulted.
    constexpr spin_delay(spin_delay const &)             noexcept = default;  //!< Defaulted.
    constexpr spin_delay(spin_delay &&)                  noexcept = default;  //!< Defaulted.
    constexpr spin_delay & operator=(spin_delay const &) noexcept = default;  //!< Defaulted.
    constexpr spin_delay & operator=(spin_delay &&)      noexcept = default;  //!< Defaulted.
    ~spin_delay()                                        noexcept = default;  //!< Defaulted.
    //!\}

    /*!\brief Delays the calling thread by either using active spinning or passive spinning.
     *
     * \details
     *
     * In the first x repetitions the thread is paused using efficient CPU instructions. After x cycles of
     * active spinning the thread will be suspended by invoking std::this_thread::yield.
     */
    void wait()
    {
        if (current <= max_repetitions)  // Start active spinning phase
        {
            for (int_fast32_t i = 0; i < current; ++i)
                pause_processor();
            current <<= 1;  // double the amount of active CPU waiting cycles.
        }
        else  // Start passive spinning phase
        {
            std::this_thread::yield();
        }
    }

private:

    //!\brief Efficient instruction to pause the CPU.
    void pause_processor()
    {
        #if defined(__SSE2__)  // AMD and Intel
            _mm_pause();
        #elif defined(__armel__) || defined(__ARMEL__) // arm, but broken? ; repeat of default case as armel also defines __arm__
            asm volatile ("nop" ::: "memory");  // default operation - does nothing => Might lead to passive spinning.
        #elif defined(__arm__) || defined(__aarch64__) // arm big endian / arm64
            __asm__ __volatile__ ("yield" ::: "memory");
        #elif defined(__ia64__)  // IA64
            __asm__ __volatile__ ("hint @pause");
        #elif defined(__powerpc__) || defined(__ppc__) || defined(__PPC__) // PowerPC
            __asm__ __volatile__ ("or 27,27,27" ::: "memory");
        #else  // everything else.
            asm volatile ("nop" ::: "memory");  // default operation - does nothing => Might lead to passive spinning.
        #endif
    }

    //!\brief The maximal number of repetitions until the thread yields.
    static constexpr int_fast32_t max_repetitions{16};
    //!\brief The current waiting phase.
    int_fast32_t                  current{1};
};

} // namespace seqan3::detail
