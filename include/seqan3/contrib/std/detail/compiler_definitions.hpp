// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan-std/blob/main/LICENSE
// -----------------------------------------------------------------------------------------------------

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
