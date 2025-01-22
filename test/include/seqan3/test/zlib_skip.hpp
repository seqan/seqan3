// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides macros for skipping tests that rely on the binary compressed output of zlib.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <gtest/gtest.h>

#include <seqan3/core/platform.hpp>

#if defined(SEQAN3_HAS_ZLIB)
#    include <zlib.h>
#endif

// Some of our tests check the binary compressed output of zlib. This is not guaranteed to be the same for all zlib
// implementations.
// This macro should be set to 1 if the zlib implementation is not the standard zlib, for example, zlib-ng.
// zlib-ng is automatically detected if `zlib.h` resolves to zlib-ng's header.
#ifndef SEQAN3_TEST_SKIP_ZLIB_DEFLATE
#    ifdef ZLIBNG_VERSION
#        define SEQAN3_TEST_SKIP_ZLIB_DEFLATE 1
#    else
#        define SEQAN3_TEST_SKIP_ZLIB_DEFLATE 0
#    endif
#endif

// Defines a GTEST_SKIP macro if the zlib implementation is not the standard zlib.
// Otherwise, it does nothing.
#ifndef SEQAN3_TEST_GTEST_SKIP_ZLIB_DEFLATE
#    if SEQAN3_TEST_SKIP_ZLIB_DEFLATE
#        define SEQAN3_TEST_GTEST_SKIP_ZLIB_DEFLATE                                                                    \
            GTEST_SKIP() << "Not testing binary compressed output for alternative zlib implementations."
#    else
#        define SEQAN3_TEST_GTEST_SKIP_ZLIB_DEFLATE
#    endif
#endif
