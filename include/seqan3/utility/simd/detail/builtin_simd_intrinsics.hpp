// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides intrinsics include for builtin simd.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

// Exclude powerpc since it may have this header and triggers a warning (-DNO_WARN_X86_INTRINSICS) which tells you that
// x86intrin.h is only there to allow porting x86_64 code to powerpc, specifically Intel intrinsics to powerpc64le.
// Since we will not support powerpc for the builtin simd backend, we will avoid including this header.
//
// See the following link for a full description of the x86intrin.h header on powerpc
// https://github.com/gcc-mirror/gcc/blob/41d6b10e96a1de98e90a7c0378437c3255814b16/gcc/config/rs6000/xmmintrin.h#L27-L55
//
// We also exclude clang with libc++ on any non-x86 platform.
// `x86intrin.h` includes `immintrin.h`, which will error:
// https://github.com/llvm/llvm-project/blob/9886788a8a500a1b429a6db64397c849b112251c/clang/lib/Headers/immintrin.h#L14
#if __has_include(<x86intrin.h>) && \
    !(defined(__powerpc__) || defined(__ppc__) || defined(__PPC__)) && \
    !(defined(__clang__) && defined(_LIBCPP_VERSION) && !defined(__i386__) && !defined(__x86_64__))
#    include <x86intrin.h> // x86 intrinsics (linux)
#endif

#if defined(_WIN32) && __has_include(<intrin.h>)
#    include <intrin.h> // x86 intrinsics (windows)
#endif

// MSVC doesn't define SSE4 macros, even if the instruction set is available (e.g. when AVX is defined)
// See https://docs.microsoft.com/en-us/cpp/preprocessor/predefined-macros?view=msvc-160
#if defined(_MSC_VER) && defined(__AVX__) && !defined(__SSE4_1__) && !defined(__SSE4_2__)
#    define __SSE4_1__ 1
#    define __SSE4_2__ 1
#endif
