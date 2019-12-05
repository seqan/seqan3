// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides intrinsics include for builtin simd.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

#if __has_include(<x86intrin.h>)
#include <x86intrin.h> // x86 intrinsics (linux)
#endif

#if __has_include(<intrin.h>)
#include <intrin.h> // x86 intrinsics (windows)
#endif

// MSVC doesn't define SSE4 macros, even if the instruction set is available (e.g. when AVX is defined)
#if defined(_MSC_VER) && defined(__AVX__) && !defined(__SSE4_1__) && !defined(__SSE4_2__)
#define __SSE4_1__ 1
#define __SSE4_2__ 1
#endif
