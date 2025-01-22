// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides platform and dependency checks.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <version>

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
#if defined(__GNUC__) && !defined(__llvm__) && !defined(__INTEL_COMPILER) && !defined(__INTEL_LLVM_COMPILER)
#    define SEQAN3_COMPILER_IS_GCC 1
#else
#    define SEQAN3_COMPILER_IS_GCC 0
#endif

#if SEQAN3_DOXYGEN_ONLY(1) 0
//!\brief This disables the warning you would get if your compiler is not known to work.
#    define SEQAN3_DISABLE_COMPILER_CHECK
#endif // SEQAN3_DOXYGEN_ONLY(1)0

// ============================================================================
//  Compiler support
// ============================================================================

#if SEQAN3_COMPILER_IS_GCC && (__GNUC__ < 12)
#    error "At least GCC 12 is needed."
#endif

// clang-format off
#if defined(__INTEL_LLVM_COMPILER) && (__INTEL_LLVM_COMPILER < 20240000)
#    error "At least Intel OneAPI 2024 is needed."
#endif
// clang-format on

#if defined(__clang__) && defined(__clang_major__) && (__clang_major__ < 17)
#    error "At least Clang 17 is needed."
#endif

// ============================================================================
//  Standard library support
// ============================================================================

#if defined(_LIBCPP_VERSION) && (_LIBCPP_VERSION < 170000)
#    error "At least libc++ 17 is required."
#endif

#if defined(_GLIBCXX_RELEASE) && (_GLIBCXX_RELEASE < 12)
#    error "At least libstdc++ 12 is needed."
#endif

// ============================================================================
//  C++ standard and features
// ============================================================================

// C++ standard [required]
#ifdef __cplusplus
#    if (__cplusplus < 202100)
#        error "SeqAn3 requires C++23, make sure that you have set -std=c++23."
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

#ifndef SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY
#    if SEQAN3_COMPILER_IS_GCC
// For checking whether workaround applies, e.g., in search_scheme_test
#        define SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY 1
// The goal is to create _Pragma("GCC diagnostic ignored \"-Wrestrict\"")
// The outer quotes are added by SEQAN3_PRAGMA, so we need SEQAN3_PRAGMA(GCC diagnostic ignored "-Wrestrict")
// SEQAN3_CONCAT_STRING(GCC diagnostic ignored, -Wrestrict) -> SEQAN3_PRAGMA(GCC diagnostic ignored "-Wrestrict")
#        define SEQAN3_CONCAT_STRING(x, y) SEQAN3_PRAGMA(x #y)
#        define SEQAN3_GCC_DIAGNOSTIC_IGNORE1(x, ...)                                                                  \
            SEQAN3_PRAGMA(GCC diagnostic push)                                                                         \
            SEQAN3_CONCAT_STRING(GCC diagnostic ignored, x)
#        define SEQAN3_GCC_DIAGNOSTIC_IGNORE2(x, y)                                                                    \
            SEQAN3_PRAGMA(GCC diagnostic push)                                                                         \
            SEQAN3_CONCAT_STRING(GCC diagnostic ignored, x)                                                            \
            SEQAN3_CONCAT_STRING(GCC diagnostic ignored, y)
// A helper that enables SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY_START to take one or two arguments
// SEQAN3_GCC_DIAGNOSTIC_IGNORE(-Wrestict, 2, 1) -> SEQAN3_GCC_DIAGNOSTIC_IGNORE1(-Wrestict, 2)
// SEQAN3_GCC_DIAGNOSTIC_IGNORE(-Wrestict, -Warray-bounds, 2, 1) -> SEQAN3_GCC_DIAGNOSTIC_IGNORE2(-Wrestict, -Warray-bounds)
#        define SEQAN3_GCC_DIAGNOSTIC_IGNORE(x, y, n, ...) SEQAN3_GCC_DIAGNOSTIC_IGNORE##n(x, y)
// SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY_START(-Wrestrict) -> SEQAN3_GCC_DIAGNOSTIC_IGNORE(-Wrestict, 2, 1)
// SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY_START(-Wrestrict, -Warray-bounds) -> SEQAN3_GCC_DIAGNOSTIC_IGNORE(-Wrestict, -Warray-bounds, 2, 1)
#        define SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY_START(x, ...) SEQAN3_GCC_DIAGNOSTIC_IGNORE(x, ##__VA_ARGS__, 2, 1)
#        define SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY_STOP SEQAN3_PRAGMA(GCC diagnostic pop)
#    else
/*!\name Workaround for bogus memcopy/memmove warnings on GCC
     * \{
     */
//!\brief Indicates whether the workaround is active. `1` for GCC, `0` for other compilers.
#        define SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY 0
/*!\brief Denotes the start of a block where diagnostics are ignored.
 * \details
 * If SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY is 0, this macro has no effect.
 * Otherwise, the macro takes one or two arguments and will expand to a preprocessor directive equivalent to:
 * ### Input
 * ```cpp
 * // The macro accepts one or two arguments.
 * SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY_START(-Wrestrict, -Warray-bounds)
 * ```
 * ### Output
 * ```cpp
 * #pragma GCC diagnostic push
 * #pragma GCC diagnostic ignored "-Wrestrict"
 * #pragma GCC diagnostic ignored "-Warray-bounds"
 * ```
 */
#        define SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY_START(...)
/*!\brief Denotes the end of a block where diagnostics are ignored.
 * \details
 * If SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY is 0, this macro has no effect.
 * Otherwise, the macro will expand to a preprocessor directive equivalent to:
 * ```cpp
 * #pragma GCC diagnostic pop
 * ```
 */
#        define SEQAN3_WORKAROUND_GCC_BOGUS_MEMCPY_STOP
//!\}
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
