// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <cinttypes>
#include <ciso646> // makes _LIBCPP_VERSION available
#include <cstddef> // makes __GLIBCXX__ available

/*!\file
 * \brief Provides platform and dependency checks.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

// macro cruft
//!\cond
#define STR_HELPER(x) #x
#define STR(x) STR_HELPER(x)
//!\endcond

// ============================================================================
//  C++ standard and features
// ============================================================================

// C++ standard [required]
#ifdef __cplusplus
    static_assert(__cplusplus >= 201703, "SeqAn3 requires C++17, make sure that you have set -std=c++17.");
#else
#   error "This is not a C++ compiler."
#endif

#if __has_include(<version>)
#   include <version>
#endif

// C++ Concepts [required]
#ifdef __cpp_concepts
#   if __cpp_concepts == 201507 // GCC and Concepts TS
#       define SEQAN3_CONCEPT concept bool
#   else
        static_assert(__cpp_concepts >= 201507, "Your compiler supports Concepts, but the support is not recent enough.");
#       define SEQAN3_CONCEPT concept
#   endif
#else
#   error "SeqAn3 requires C++ Concepts, either vie the TS (flag: -fconcepts) or via C++20 (flag: -std=c++2a / -std=c++20)."
#endif

// filesystem [required]
#if !__has_include(<filesystem>)
#   if !__has_include(<experimental/filesystem>)
#      error SeqAn3 requires C++17 filesystem support, but it was not found.
#   endif
#endif

// ============================================================================
//  Dependencies
// ============================================================================

// SeqAn [required]
#if !__has_include(<seqan3/version.hpp>)
#   error SeqAn3 include directory not set correctly. Forgot to add -I ${INSTALLDIR}/include to your CXXFLAGS?
#endif

// Ranges [required]
#if __has_include(<range/v3/version.hpp>)
#   define RANGE_V3_MINVERSION 1000
#   define RANGE_V3_MAXVERSION 1099
// TODO the following doesn't actually show the current version, only its formula. How'd you do it?
#   define MSG "Your version: " STR(RANGE_V3_VERSION) \
                "; minimum version: " STR(RANGE_V3_MINVERSION) \
                "; expected maximum version: " STR(RANGE_V3_MAXVERSION)
#   include <range/v3/version.hpp>
#   if RANGE_V3_VERSION < RANGE_V3_MINVERSION
#       error Your range-v3 library is too old.
#       pragma message(MSG)
#   elif RANGE_V3_VERSION > RANGE_V3_MAXVERSION
#       pragma GCC warning "Your range-v3 library is possibly too new. Some features might not work correctly."
#       pragma message(MSG)
#   endif
#   undef MSG
#else
#   error The range-v3 library was not included correctly. Forgot to add -I ${INSTALLDIR}/include to your CXXFLAGS?
#endif

// SDSL [required]
#if __has_include(<sdsl/version.hpp>)
#   include <sdsl/version.hpp>
    static_assert(sdsl::sdsl_version_major == 3, "Only version 3 of the SDSL is supported by SeqAn3.");
#else
#   error The sdsl library was not included correctly. Forgot to add -I ${INSTALLDIR}/include to your CXXFLAGS?
#endif

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

     //! \brief Macro for Cereal's serialize function.
#   define CEREAL_SERIALIZE_FUNCTION_NAME serialize
     //! \brief Macro for Cereal's load function.
#   define CEREAL_LOAD_FUNCTION_NAME load
     //! \brief Macro for Cereal's save function.
#   define CEREAL_SAVE_FUNCTION_NAME save
     //! \brief Macro for Cereal's load_minimal function.
#   define CEREAL_LOAD_MINIMAL_FUNCTION_NAME load_minimal
     //! \brief Macro for Cereal's save_minimal function.
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

// ============================================================================
//  Documentation
// ============================================================================

// Doxygen related
// this macro is a NO-OP unless doxygen parses it, in which case it resolves to the argument
#ifndef SEQAN3_DOXYGEN_ONLY
#   define SEQAN3_DOXYGEN_ONLY(x)
#endif

// ============================================================================
//  Deprecation Messages
// ============================================================================

//!\brief Deprecation message for SeqAn 3.1.0 release.
#if !defined(SEQAN3_DEPRECATED_310)
#   define SEQAN3_DEPRECATED_310 [[deprecated("This will be removed in SeqAn-3.1.0; please see the documentation.")]]
#endif

// ============================================================================
//  Workarounds
// ============================================================================

#ifndef SEQAN3_WORKAROUND_VIEW_PERFORMANCE
//!\brief Performance of views, especially filter and join is currently bad, especially in I/O.
#   define SEQAN3_WORKAROUND_VIEW_PERFORMANCE 1
#endif

//!\brief See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=87113
#ifndef SEQAN3_WORKAROUND_GCC_87113
#   if defined(__GNUC_MINOR__) && ((__GNUC__ == 7) || ((__GNUC__ == 8) && (__GNUC_MINOR__ < 3)))
#       define SEQAN3_WORKAROUND_GCC_87113 1
#   else
#       define SEQAN3_WORKAROUND_GCC_87113 0
#   endif
#endif

//!\brief See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90897
#ifndef SEQAN3_WORKAROUND_GCC_90897
#   if defined(__GNUC__) && (__GNUC__ == 8)
#       define SEQAN3_WORKAROUND_GCC_90897 1
#   else
#       define SEQAN3_WORKAROUND_GCC_90897 0
#   endif
#endif

//!\brief Various concept problems only present in GCC7 and GCC8.
#ifndef SEQAN3_WORKAROUND_GCC7_AND_8_CONCEPT_ISSUES
#   if defined(__GNUC__) && ((__GNUC__ == 7) || (__GNUC__ == 8))
#       define SEQAN3_WORKAROUND_GCC7_AND_8_CONCEPT_ISSUES 1
#   else
#       define SEQAN3_WORKAROUND_GCC7_AND_8_CONCEPT_ISSUES 0
#   endif
#endif

// ============================================================================
//  Backmatter
// ============================================================================

// macro cruft undefine
#undef STR
#undef STR_HELPER
