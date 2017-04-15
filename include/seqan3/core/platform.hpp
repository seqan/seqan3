// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

#pragma once

#include <cinttypes>
#include <ciso646> // makes _LIBCPP_VERSION available
#include <cstddef> // makes __GLIBCXX__ available

/*!\file platform.hpp
 * \brief Contains platform and dependency checks.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \ingroup core
 */

// macro cruft
//!\cond
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
//!\endcond

// C++ standard
#ifdef __cplusplus
    static_assert(__cplusplus >= 201500, "SeqAn3 requires C++17, make sure that you have set -std=c++17.");
#else
#   error "This is not a C++ compiler."
#endif

// Concepts TS
#ifdef __cpp_concepts
    static_assert(__cpp_concepts >= 201507, "Your compiler supports Concepts, but the support is not recent enough.");
#else
#   error "SeqAn3 requires the Concepts TS, make sure that you have set -fconcepts (not all compilers support this)."
#endif

// SeqAn
#if !__has_include(<seqan3/version.hpp>)
#   error SeqAn3 include directory not set correctly. Forgot to add -I ${INSTALLDIR}/include to your CXXFLAGS?
#endif

// Ranges
#if __has_include(<range/v3/version.hpp>)
#   define RANGE_V3_MINVERSION 200
#   define RANGE_V3_MAXVERSION 299
// TODO the following doesn't actually show the current version, only its formula. How'd you do it?
#   define MSG "Your version: " STR(RANGE_V3_VERSION) \
                "; minimum version: " STR(RANGE_V3_MINVERSION) \
                "; expected maximum version: " STR(RANGE_V3_MAXVERSION)
#   include <range/v3/version.hpp>
#   if RANGE_V3_VERSION < RANGE_V3_MINVERSION
#       error Your range-v3 library is too old.
#       pragma message(MSG)
#   elif RANGE_V3_VERSION > RANGE_V3_MAXVERSION
#       pragma GCC warning "Your range-v3 library is possibly tot new. Some features might not work correctly."
#       pragma message(MSG)
#   endif
#   undef MSG
#else
#   error The range-v3 library was not included correctly. Forgot to add -I ${INSTALLDIR}/include to your CXXFLAGS?
#endif

// SDSL
// TODO (doesn't have a version.hpp, yet)

#undef STR
#undef STR_HELPER
