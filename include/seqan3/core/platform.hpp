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

/*!\file
 * \brief Contains platform and dependency checks.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

// macro cruft
//!\cond
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
//!\endcond

// C++ standard [required]
#ifdef __cplusplus
    static_assert(__cplusplus >= 201500, "SeqAn3 requires C++17, make sure that you have set -std=c++17.");
#else
#   error "This is not a C++ compiler."
#endif

// Concepts TS [required]
#ifdef __cpp_concepts
    static_assert(__cpp_concepts >= 201507, "Your compiler supports Concepts, but the support is not recent enough.");
#else
#   error "SeqAn3 requires the Concepts TS, make sure that you have set -fconcepts (not all compilers support this)."
#endif

// SeqAn [required]
#if !__has_include(<seqan3/version.hpp>)
#   error SeqAn3 include directory not set correctly. Forgot to add -I ${INSTALLDIR}/include to your CXXFLAGS?
#endif

// Ranges [required]
#if __has_include(<range/v3/version.hpp>)
#   define RANGE_V3_MINVERSION 300
#   define RANGE_V3_MAXVERSION 399
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

// SDSL [required]
// TODO (doesn't have a version.hpp, yet)

// Cereal [optional]
/*!\def SEQAN3_WITH_CEREAL
 * \brief Whether CEREAL support is available or not.
 * \ingroup core
 */
#ifndef SEQAN3_WITH_CEREAL
#   if __has_include(<cereal/cereal.hpp>)
#       define SEQAN3_WITH_CEREAL 1
#   else
#       define SEQAN3_WITH_CEREAL 0
#   endif
#elif SEQAN3_WITH_CEREAL != 0
#   if ! __has_include(<cereal/cereal.hpp>)
#       error Cereal was marked as required, but not found!
#   endif
#endif

#if !SEQAN3_WITH_CEREAL
    /*!\cond DEV
     * \name Cereal function macros
     * \ingroup core
     * \brief These can be changed by apps so we used the macros instead of the values internally.
     * \{
     */
#   define CEREAL_SERIALIZE_FUNCTION_NAME serialize
#   define CEREAL_LOAD_FUNCTION_NAME load
#   define CEREAL_SAVE_FUNCTION_NAME save
#   define CEREAL_LOAD_MINIMAL_FUNCTION_NAME load_minimal
#   define CEREAL_SAVE_MINIMAL_FUNCTION_NAME save_minimal
    /*!\}
     * \endcond
     */
#endif

// Lemon [optional]
/*!\def SEQAN3_WITH_LEMON
 * \brief Whether Lemon support is available or not.
 * \ingroup core
 */
#ifndef SEQAN3_WITH_LEMON
#   if __has_include(<lemon/config.h>)
#       define SEQAN3_WITH_LEMON 1
#   else
#       define SEQAN3_WITH_LEMON 0
#   endif
#elif SEQAN3_WITH_LEMON != 0
#   if !__has_include(<lemon/config.h>)
#       error Lemon was marked as required, but not found!
#   endif
#endif
#if SEQAN3_WITH_LEMON == 1
#   define LEMON_HAVE_LONG_LONG 1
#   define LEMON_CXX11 1
#   if defined(__unix__) || defined(__APPLE__)
#       define LEMON_USE_PTHREAD 1
#       define LEMON_USE_WIN32_THREADS 0
#       define LEMON_WIN32 0
#   else
#       define LEMON_USE_PTHREAD 0
#       define LEMON_USE_WIN32_THREADS 1
#       define LEMON_WIN32 1
#   endif
#endif

// TODO (doesn't have a version.hpp, yet)

#undef STR
#undef STR_HELPER
