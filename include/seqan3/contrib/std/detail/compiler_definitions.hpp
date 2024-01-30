// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides compiler definitions for seqan::stl.
 */

// File might be included from multiple libraries.
#ifndef SEQAN_STD_DETAIL_COMPILER_DEFINITIONS
#define SEQAN_STD_DETAIL_COMPILER_DEFINITIONS

#if defined(__GNUC__) && !defined(__llvm__) && !defined(__INTEL_COMPILER)
#    define SEQAN_STD_COMPILER_IS_GCC 1
#else
#    define SEQAN_STD_COMPILER_IS_GCC 0
#endif

// Bug: Nested classes are not friends of outer classes.
// Fixed with GCC13.
#if SEQAN_STD_COMPILER_IS_GCC && (__GNUC__ > 12)
#    define SEQAN_STD_NESTED_VISIBILITY private:
#else
#    define SEQAN_STD_NESTED_VISIBILITY public:
#endif

#endif // SEQAN_STD_DETAIL_COMPILER_DEFINITIONS
