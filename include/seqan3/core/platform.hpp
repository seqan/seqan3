// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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
//  Compiler support general
// ============================================================================

/*!\def SEQAN3_COMPILER_IS_GCC
 * \brief Whether the current compiler is GCC.
 * \ingroup core
 * \details
 * __GNUC__ is also used to indicate the support for GNU compiler extensions. To detect the presence of the GCC
 * compiler, one has to rule out other compilers.
 *
 * \sa https://sourceforge.net/p/predef/wiki/Compilers
 */
#if defined(__GNUC__) && !defined(__llvm__) && !defined(__INTEL_COMPILER)
#    define SEQAN3_COMPILER_IS_GCC 1
#else
#    define SEQAN3_COMPILER_IS_GCC 0
#endif

#if SEQAN3_DOXYGEN_ONLY(1) 0
//!\brief This disables the warning you would get if your compiler is not known to work.
#    define SEQAN3_DISABLE_COMPILER_CHECK
#endif // SEQAN3_DOXYGEN_ONLY(1)0

// ============================================================================
//  Compiler support GCC
// ============================================================================

#if SEQAN3_COMPILER_IS_GCC
#    if (__GNUC__ < 11)
#        error                                                                                                         \
            "SeqAn 3.1.x is the last version that supports GCC 7, 8, and 9. SeqAn 3.2.x is the latest version that support GCC 10. Please upgrade your compiler or use 3.1.x./3.2.x."
#    endif // (__GNUC__ < 11)

#    if (__GNUC__ == 11 && __GNUC_MINOR__ <= 3)
#        pragma GCC warning "Be aware that GCC < 11.4 might have bugs that cause SeqAn3 fail to compile."
#    endif // (__GNUC__ == 11 && __GNUC_MINOR__ <= 2)

#    if (__GNUC__ == 12 && __GNUC_MINOR__ <= 2)
#        pragma GCC warning "Be aware that GCC < 12.3 might have bugs that cause SeqAn3 fail to compile."
#    endif // (__GNUC__ == 12 && __GNUC_MINOR__ <= 1)

#    if SEQAN3_DOXYGEN_ONLY(1) 0
//!\brief This disables the warning you would get if your compiler is newer than the latest supported version.
#        define SEQAN3_DISABLE_NEWER_COMPILER_DIAGNOSTIC
#    endif // SEQAN3_DOXYGEN_ONLY(1)0

#    ifndef SEQAN3_DISABLE_NEWER_COMPILER_DIAGNOSTIC
#        if (__GNUC__ > 14)
#            pragma message                                                                                            \
                "Your compiler is newer than the latest supported compiler of this SeqAn version (gcc-13). It might be that SeqAn does not compile due to this. You can disable this warning by setting -DSEQAN3_DISABLE_NEWER_COMPILER_DIAGNOSTIC."
#        endif // (__GNUC__ > 13)
#    endif     // SEQAN3_DISABLE_NEWER_COMPILER_DIAGNOSTIC

// ============================================================================
//  Compiler support Clang
// ============================================================================

#elif defined(__clang__)
#    if __clang_major__ < 17
#        error "Only Clang >= 17 is supported."
#    endif

// ============================================================================
//  Compiler support other
// ============================================================================

#elif !defined(SEQAN3_DISABLE_COMPILER_CHECK)
#    error "Your compiler is not supported. You can disable this error by setting -DSEQAN3_DISABLE_COMPILER_CHECK."
#endif // SEQAN3_COMPILER_IS_GCC

// ============================================================================
//  C++ standard and features
// ============================================================================

#if __has_include(<version>)
#    include <version>
#endif

// C++ standard [required]
#ifdef __cplusplus
#    if (__cplusplus < 202002L)
#        error "SeqAn3 requires C++20, make sure that you have set -std=c++20."
#    endif
#else
#    error "This is not a C++ compiler."
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

// ============================================================================
//  Deprecation Messages
// ============================================================================

//!\brief _Pragma requires a string-literal and # makes it a string
#ifndef SEQAN3_PRAGMA
#    define SEQAN3_PRAGMA(non_string_literal) _Pragma(#non_string_literal)
#endif

//!\brief Deprecation message for deprecated header.
#ifndef SEQAN3_DEPRECATED_HEADER
#    ifndef SEQAN3_DISABLE_DEPRECATED_WARNINGS
#        define SEQAN3_DEPRECATED_HEADER(message) SEQAN3_PRAGMA(GCC warning message)
#    else
#        define SEQAN3_DEPRECATED_HEADER(message) /**/
#    endif
#endif

//!\brief Deprecation message for SeqAn 3.4.0 release.
#ifndef SEQAN3_REMOVE_DEPRECATED_340
#    ifndef SEQAN3_DEPRECATED_340
#        ifndef SEQAN3_DISABLE_DEPRECATED_WARNINGS
#            define SEQAN3_DEPRECATED_340                                                                              \
                [[deprecated("This will be removed in SeqAn-3.4.0; please see the documentation.")]]
#        else
#            define SEQAN3_DEPRECATED_340 /**/
#        endif
#    endif
#endif

// ============================================================================
//  Workarounds
// ============================================================================

#ifndef SEQAN3_WORKAROUND_VIEW_PERFORMANCE
//!\brief Performance of views, especially filter and join is currently bad, especially in I/O.
#    define SEQAN3_WORKAROUND_VIEW_PERFORMANCE 1
#endif

//!\brief A view does not need to be default constructible. This change is first implemented in gcc12.
#ifndef SEQAN3_WORKAROUND_DEFAULT_CONSTRUCTIBLE_VIEW
#    if SEQAN3_COMPILER_IS_GCC && (__GNUC__ < 12)
#        define SEQAN3_WORKAROUND_DEFAULT_CONSTRUCTIBLE_VIEW 1
#    else
#        define SEQAN3_WORKAROUND_DEFAULT_CONSTRUCTIBLE_VIEW 0
#    endif
#endif

//!\brief See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=100139
//!       std::views::{take, drop} do not type-erase. This is a defect within the standard lib (fixed in gcc12).
#ifndef SEQAN3_WORKAROUND_GCC_100139
#    if SEQAN3_COMPILER_IS_GCC && (__GNUC__ < 12)
#        define SEQAN3_WORKAROUND_GCC_100139 1
#    else
#        define SEQAN3_WORKAROUND_GCC_100139 0
#    endif
#endif

/*!\brief Workaround bogus memcpy errors in GCC > 12. (Wrestrict and Wstringop-overflow)
 * \see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=105545
 */
#ifndef SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY
#    if SEQAN3_COMPILER_IS_GCC && (__GNUC__ >= 12)
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

/*!\brief Our char literals returning std::vector should be constexpr if constexpr std::vector is supported.
 *
 * The _GLIBCXX_DEBUG statement is a workaround for a libstdc++ bug
 * \see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=104748
 * \see https://github.com/seqan/seqan3/issues/3221
 * \see https://godbolt.org/z/159n8xrdo
 **/
#if __cpp_lib_constexpr_vector >= 201907L && (defined(_LIBCPP_VERSION) || !defined(_GLIBCXX_DEBUG))
#    define SEQAN3_WORKAROUND_LITERAL constexpr
#else
#    define SEQAN3_WORKAROUND_LITERAL inline
#endif

#if SEQAN3_DOXYGEN_ONLY(1) 0
//!\brief This disables the warning you would get if -D_GLIBCXX_USE_CXX11_ABI=0 is set.
#    define SEQAN3_DISABLE_LEGACY_STD_DIAGNOSTIC
#endif // SEQAN3_DOXYGEN_ONLY(1)0

#if defined(_GLIBCXX_USE_CXX11_ABI) && _GLIBCXX_USE_CXX11_ABI == 0
#    ifndef SEQAN3_DISABLE_LEGACY_STD_DIAGNOSTIC
#        pragma message                                                                                                \
            "We do not actively support compiler that have -D_GLIBCXX_USE_CXX11_ABI=0 set, and it might be that SeqAn does not compile due to this. It is known that all compiler of CentOS 7 / RHEL 7 set this flag by default (and that it cannot be overridden!). Note that these versions of the OSes are community-supported (see https://docs.seqan.de/seqan3/main_user/about_api.html#platform_stability for more details). You can disable this warning by setting -DSEQAN3_DISABLE_LEGACY_STD_DIAGNOSTIC."
#    endif // SEQAN3_DISABLE_LEGACY_STD_DIAGNOSTIC
#endif     // _GLIBCXX_USE_CXX11_ABI == 0

// ============================================================================
//  Backmatter
// ============================================================================

// macro cruft undefine
#undef SEQAN3_STR
#undef SEQAN3_STR_HELPER
