// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::spin_delay.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <thread>
#include <xmmintrin.h>

#include <seqan3/core/platform.hpp>

namespace seqan3::contrib
{
//!\cond
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

    void wait()
    {
        if (current <= max_repetitions)  // Start active spinning phase
        {
            for (int_fast32_t i = 0; i < current; ++i)
                pause_processor();
            current <<= 1;  // double the amount of active CPU waiting.
        }
        else  // Start passive spinning phase
        {
            std::this_thread::yield();
        }
    }

private:

    void pause_processor()
    {
        #if defined(__SSE2__)  // AMD and Intel
            _mm_pause();
        #elif defined(__ia64__)  // IA64
            __asm__ __volatile__ ("hint @pause");
        #else  // everything else.
            asm volatile ("nop" ::: "memory");  // default operation - does nothing => Might lead to passive spinning.
        #endif
    }

    static constexpr int_fast32_t max_repetitions{16};
    int_fast32_t                  current{1};
};
//!\endcond
} // namespace seqan3::contrib
