// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides platform and dependency checks.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cinttypes>
#include <ciso646> // makes _LIBCPP_VERSION available
#include <cstddef> // makes __GLIBCXX__ available

// macro cruft
//!\cond
#define SEQAN3_STR_HELPER(x) #x
#define SEQAN3_STR(x) SEQAN3_STR_HELPER(x)
//!\endcond

// ============================================================================
//  Documentation
// ============================================================================

// Doxygen related
// this macro is a NO-OP unless doxygen parses it, in which case it resolves to the argument
#ifndef SEQAN3_DOXYGEN_ONLY
#    define SEQAN3_DOXYGEN_ONLY(x)
#endif

// ============================================================================
//  Compiler support
// ============================================================================

#if defined(__GNUC__) && (__GNUC__ < 10)
#    error "SeqAn 3.1.x is the last version that supports GCC 7, 8, and 9. Please upgrade your compiler or use 3.1.x."
#endif // defined(__GNUC__) && (__GNUC__ < 10)

#if SEQAN3_DOXYGEN_ONLY(1) 0
//!\brief This disables the warning you would get if your compiler is newer than the latest supported version.
#    define SEQAN3_DISABLE_NEWER_COMPILER_DIAGNOSTIC
#endif // SEQAN3_DOXYGEN_ONLY(1)0

#ifndef SEQAN3_DISABLE_NEWER_COMPILER_DIAGNOSTIC
#    if defined(__GNUC__) && (__GNUC__ > 12)
#        pragma message                                                                                                \
            "Your compiler is newer than the latest supported compiler of this SeqAn version (gcc-12). It might be that SeqAn does not compile due to this. You can disable this warning by setting -DSEQAN3_DISABLE_NEWER_COMPILER_DIAGNOSTIC."
#    endif // defined(__GNUC__) && (__GNUC__ > 12)
#endif     // SEQAN3_DISABLE_NEWER_COMPILER_DIAGNOSTIC

// ============================================================================
//  C++ standard and features
// ============================================================================

// C++ standard [required]
// Note: gcc10 -std=c++20 still defines __cplusplus=201709
#ifdef __cplusplus
#    if (__cplusplus < 201709)
#        error "SeqAn3 requires C++20, make sure that you have set -std=c++20."
#    endif
#else
#    error "This is not a C++ compiler."
#endif

#if __has_include(<version>)
#    include <version>
#endif

// ============================================================================
//  Dependencies
// ============================================================================

// SeqAn [required]
#if __has_include(<seqan3/version.hpp>)
#    include <seqan3/version.hpp>
#else
#    error SeqAn3 include directory not set correctly. Forgot to add -I ${INSTALLDIR}/include to your CXXFLAGS?
#endif

// SDSL [required]
#if __has_include(<sdsl/version.hpp>)
#    include <sdsl/version.hpp>
static_assert(sdsl::sdsl_version_major == 3, "Only version 3 of the SDSL is supported by SeqAn3.");
#else
#    error The sdsl library was not included correctly. Forgot to add -I ${INSTALLDIR}/include to your CXXFLAGS?
#endif

// Cereal [optional]
/*!\def SEQAN3_WITH_CEREAL
 * \brief Whether CEREAL support is available or not.
 * \ingroup core
 */
#ifndef SEQAN3_WITH_CEREAL
#    if __has_include(<cereal/cereal.hpp>)
#        define SEQAN3_WITH_CEREAL 1
#    else
#        define SEQAN3_WITH_CEREAL 0
#    endif
#elif SEQAN3_WITH_CEREAL != 0
#    if !__has_include(<cereal/cereal.hpp>)
#        error Cereal was marked as required, but not found!
#    endif
#endif

//!\cond DEV
#if !SEQAN3_WITH_CEREAL
/*!\name Cereal function macros
 * \ingroup core
 * \brief These can be changed by apps so we used the macros instead of the values internally.
 * \{
 */
#    define CEREAL_SERIALIZE_FUNCTION_NAME serialize       //!< Macro for Cereal's serialize function.
#    define CEREAL_LOAD_FUNCTION_NAME load                 //!< Macro for Cereal's load function.
#    define CEREAL_SAVE_FUNCTION_NAME save                 //!< Macro for Cereal's save function.
#    define CEREAL_LOAD_MINIMAL_FUNCTION_NAME load_minimal //!< Macro for Cereal's load_minimal function.
#    define CEREAL_SAVE_MINIMAL_FUNCTION_NAME save_minimal //!< Macro for Cereal's save_minimal function.
//!\}
#endif
//!\endcond

// Lemon [optional]
/*!\def SEQAN3_WITH_LEMON
 * \brief Whether Lemon support is available or not.
 * \ingroup core
 */
#ifndef SEQAN3_WITH_LEMON
#    if __has_include(<lemon/config.h>)
#        define SEQAN3_WITH_LEMON 1
#    else
#        define SEQAN3_WITH_LEMON 0
#    endif
#elif SEQAN3_WITH_LEMON != 0
#    if !__has_include(<lemon/config.h>)
#        error Lemon was marked as required, but not found!
#    endif
#endif
#if SEQAN3_WITH_LEMON == 1
#    define LEMON_HAVE_LONG_LONG 1
#    define LEMON_CXX11 1
#    if defined(__unix__) || defined(__APPLE__)
#        define LEMON_USE_PTHREAD 1
#        define LEMON_USE_WIN32_THREADS 0
#        define LEMON_WIN32 0
#    else
#        define LEMON_USE_PTHREAD 0
#        define LEMON_USE_WIN32_THREADS 1
#        define LEMON_WIN32 1
#    endif
#endif

// TODO (doesn't have a version.hpp, yet)

// ============================================================================
//  Deprecation Messages
// ============================================================================

//!\brief _Pragma requires a string-literal and # makes it a string
#ifndef SEQAN3_PRAGMA
#    define SEQAN3_PRAGMA(non_string_literal) _Pragma(#    non_string_literal)
#endif

//!\brief Deprecation message for deprecated header.
#ifndef SEQAN3_DEPRECATED_HEADER
#    ifndef SEQAN3_DISABLE_DEPRECATED_WARNINGS
#        define SEQAN3_DEPRECATED_HEADER(message) SEQAN3_PRAGMA(GCC warning message)
#    else
#        define SEQAN3_DEPRECATED_HEADER(message) /**/
#    endif
#endif

//!\brief Deprecation message for SeqAn 3.3.0 release.
#ifndef SEQAN3_REMOVE_DEPRECATED_330
#    ifndef SEQAN3_DEPRECATED_330
#        ifndef SEQAN3_DISABLE_DEPRECATED_WARNINGS
#            define SEQAN3_DEPRECATED_330                                                                              \
                [[deprecated("This will be removed in SeqAn-3.3.0; please see the documentation.")]]
#        else
#            define SEQAN3_DEPRECATED_330 /**/
#        endif
#    endif
#endif

// ============================================================================
//  Workarounds
// ============================================================================

/*!\brief Warn about gcc 10.0 and gcc 10.1
 * Known Issues:
 * * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=93983
 * * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95371
 * * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=95578
 * * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=99318
 * * https://gcc.gnu.org/bugzilla/show_bug.cgi?id=93467
 */
#if defined(__GNUC__) && (__GNUC__ == 10 && __GNUC_MINOR__ <= 2)
#    pragma GCC warning                                                                                                \
        "Be aware that gcc 10.0, 10.1 and 10.2 are known to have several bugs that might make SeqAn3 fail to compile. Please use gcc >= 10.3."
#endif // defined(__GNUC__) && (__GNUC__ == 10 && __GNUC_MINOR__ <= 1)

#ifndef SEQAN3_WORKAROUND_VIEW_PERFORMANCE
//!\brief Performance of views, especially filter and join is currently bad, especially in I/O.
#    define SEQAN3_WORKAROUND_VIEW_PERFORMANCE 1
#endif

//!\brief See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=96070 and https://github.com/seqan/product_backlog/issues/151
#ifndef SEQAN3_WORKAROUND_GCC_96070 // fixed since gcc10.4
#    if defined(__GNUC__) && (__GNUC__ == 10 && __GNUC_MINOR__ < 4)
#        define SEQAN3_WORKAROUND_GCC_96070 1
#    else
#        define SEQAN3_WORKAROUND_GCC_96070 0
#    endif
#endif

//!\brief Workaround to access the static id of the configuration elements inside of the concept definition
//!       (fixed in gcc11).
#ifndef SEQAN3_WORKAROUND_GCC_PIPEABLE_CONFIG_CONCEPT
#    if defined(__GNUC__) && (__GNUC__ < 11)
#        define SEQAN3_WORKAROUND_GCC_PIPEABLE_CONFIG_CONCEPT 1
#    else
#        define SEQAN3_WORKAROUND_GCC_PIPEABLE_CONFIG_CONCEPT 0
#    endif
#endif

//!\brief A view does not need to be default constructible. This change is first implemented in gcc12.
#ifndef SEQAN3_WORKAROUND_DEFAULT_CONSTRUCTIBLE_VIEW
#    if defined(__GNUC__) && (__GNUC__ < 12)
#        define SEQAN3_WORKAROUND_DEFAULT_CONSTRUCTIBLE_VIEW 1
#    else
#        define SEQAN3_WORKAROUND_DEFAULT_CONSTRUCTIBLE_VIEW 0
#    endif
#endif

//!\brief See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=100139
//!       std::views::{take, drop} do not type-erase. This is a defect within the standard lib (fixed in gcc12).
#ifndef SEQAN3_WORKAROUND_GCC_100139
#    if defined(__GNUC__) && (__GNUC__ < 12)
#        define SEQAN3_WORKAROUND_GCC_100139 1
#    else
#        define SEQAN3_WORKAROUND_GCC_100139 0
#    endif
#endif

/*!\brief Workaround bogus memcpy errors in GCC 12.1. (Wrestrict and Wstringop-overflow)
 * \see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=105545
 */
#ifndef SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY
#    if defined(__GNUC__) && (__GNUC__ == 12 && __GNUC_MINOR__ < 2)
#        define SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY 1
#    else
#        define SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY 0
#    endif
#endif

/*!\brief This is needed to support CentOS 7 or RHEL 7; Newer CentOS's include a more modern default-gcc version making
 *        this macro obsolete.
 *
 * In GCC 5 there was a bigger ABI change and modern systems compile with dual ABI, but some enterprise systems (those
 * where gcc 4 is the standard compiler) don't support dual ABI. This has the effect that even community builds of gcc
 * are build with --disable-libstdcxx-dual-abi. Only building the compiler yourself would solve this problem.
 *
 * \see https://github.com/seqan/seqan3/issues/2244
 * \see https://gcc.gnu.org/onlinedocs/libstdc++/manual/using_dual_abi.html
 */
#ifndef SEQAN3_WORKAROUND_GCC_NO_CXX11_ABI
#    if defined(_GLIBCXX_USE_CXX11_ABI) && _GLIBCXX_USE_CXX11_ABI == 0
#        define SEQAN3_WORKAROUND_GCC_NO_CXX11_ABI 1
#    else
#        define SEQAN3_WORKAROUND_GCC_NO_CXX11_ABI 0
#    endif
#endif

#if SEQAN3_DOXYGEN_ONLY(1) 0
//!\brief This disables the warning you would get if -D_GLIBCXX_USE_CXX11_ABI=0 is set.
#    define SEQAN3_DISABLE_LEGACY_STD_DIAGNOSTIC
#endif // SEQAN3_DOXYGEN_ONLY(1)0

#if defined(_GLIBCXX_USE_CXX11_ABI) && _GLIBCXX_USE_CXX11_ABI == 0
#    ifndef SEQAN3_DISABLE_LEGACY_STD_DIAGNOSTIC
#        pragma message                                                                                                \
            "We do not actively support compiler that have -D_GLIBCXX_USE_CXX11_ABI=0 set, and it might be that SeqAn does not compile due to this. It is known that all compiler of CentOS 7 / RHEL 7 set this flag by default (and that it cannot be overridden!). Note that these versions of the OSes are community-supported (see https://docs.seqan.de/seqan/3-master-user/about_api.html#platform_stability for more details). You can disable this warning by setting -DSEQAN3_DISABLE_LEGACY_STD_DIAGNOSTIC."
#    endif // SEQAN3_DISABLE_LEGACY_STD_DIAGNOSTIC
#endif     // _GLIBCXX_USE_CXX11_ABI == 0

// ============================================================================
//  Backmatter
// ============================================================================

// macro cruft undefine
#undef SEQAN3_STR
#undef SEQAN3_STR_HELPER

// ============================================================================
// API - Adds deprecations for seqan3/std/iterator. Will be removed in 3.3.0.
// This is the only file that is guaranteed to be included.
// Deprecations cannot be in previous header because changing the include to
// <iterator> will not include the cpp20 definitions.
// ============================================================================

//!\cond
namespace std
{

template <class Container>
class back_insert_iterator;

template <class Container>
constexpr back_insert_iterator<Container> back_inserter(Container & c);

template <class CharT>
class char_traits;

template <class T, class CharT, class Traits>
class ostream_iterator;

template <class CharT, class Traits>
class ostreambuf_iterator;

} // namespace std

namespace std::cpp20
{

// Extra include guard is needed for header tests. Prevents redefinition.
#ifndef SEQAN3_CPP20_ODR
#    define SEQAN3_CPP20_ODR 1
template <class Container>
SEQAN3_DEPRECATED_330 inline constexpr std::back_insert_iterator<Container> back_inserter(Container & c)
{
    return std::back_inserter(c);
}
#endif

template <class T, class CharT = char, class Traits = std::char_traits<CharT>>
using ostream_iterator SEQAN3_DEPRECATED_330 = std::ostream_iterator<T, CharT, Traits>;

template <class CharT, class Traits = std::char_traits<CharT>>
using ostreambuf_iterator SEQAN3_DEPRECATED_330 = std::ostreambuf_iterator<CharT, Traits>;

} // namespace std::cpp20
//!\endcond
