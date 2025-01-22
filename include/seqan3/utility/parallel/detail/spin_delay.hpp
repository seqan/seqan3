// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::spin_delay.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <thread>

#include <seqan3/core/platform.hpp>

//!\cond
#ifndef SEQAN3_HAS_MM_PAUSE
#    if defined(__SSE2__) && __has_include(<xmmintrin.h>)
#        include <xmmintrin.h> // _mm_pause()
#        define SEQAN3_HAS_MM_PAUSE 1
#    else
#        define SEQAN3_HAS_MM_PAUSE 0
#    endif // defined(__SSE2__) && __has_include(<xmmintrin.h>)
#endif     // SEQAN3_HAS_MM_PAUSE
//!\endcond

namespace seqan3::detail
{
/*!\brief A delay for threads waiting for a shared resource.
 * \ingroup utility_parallel
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
    constexpr spin_delay() noexcept = default;                               //!< Defaulted.
    constexpr spin_delay(spin_delay const &) noexcept = default;             //!< Defaulted.
    constexpr spin_delay(spin_delay &&) noexcept = default;                  //!< Defaulted.
    constexpr spin_delay & operator=(spin_delay const &) noexcept = default; //!< Defaulted.
    constexpr spin_delay & operator=(spin_delay &&) noexcept = default;      //!< Defaulted.
    ~spin_delay() noexcept = default;                                        //!< Defaulted.

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
        if (current <= max_repetitions) // Start active spinning phase
        {
            for (int_fast32_t i = 0; i < current; ++i)
                pause_processor();
            current <<= 1; // double the amount of active CPU waiting cycles.
        }
        else // Start passive spinning phase
        {
            std::this_thread::yield();
        }
    }

private:
    //!\brief Efficient instruction to pause the CPU.
    void pause_processor()
    {
#if SEQAN3_HAS_MM_PAUSE // AMD and Intel
        _mm_pause();
#elif defined(__armel__) || defined(__ARMEL__) // ARM, but broken? repeat of default case as ARMEL also defines __arm__
        asm volatile("nop" ::: "memory"); // default operation - does nothing => Might lead to passive spinning.
#elif defined(__arm__) || defined(__aarch64__) // ARM big endian / ARM64
        __asm__ __volatile__("yield" ::: "memory");
#elif defined(__ia64__)                        // IA64
        __asm__ __volatile__("hint @pause");
#elif defined(__powerpc__) || defined(__ppc__) || defined(__PPC__) || defined(__ppc64__) // PowerPC
#    if defined(__APPLE__)
        __asm__ volatile("or r27,r27,r27" ::: "memory");
#    else
        __asm__ __volatile__("or 27,27,27" ::: "memory");
#    endif
#else // everything else
        asm volatile("nop" ::: "memory"); // default operation - does nothing => Might lead to passive spinning.
#endif
    }

    //!\brief The maximal number of repetitions until the thread yields.
    static constexpr int_fast32_t max_repetitions{16};
    //!\brief The current waiting phase.
    int_fast32_t current{1};
};

} // namespace seqan3::detail
