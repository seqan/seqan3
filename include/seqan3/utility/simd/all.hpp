// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Meta-header for the \link utility_simd Utility / SIMD submodule \endlink.
 */

/*!\defgroup utility_simd SIMD
 * \brief The simd module contains a unified interface to provide simd types and functions used in seqan3.
 * \ingroup utility
 * \see utility
 * \see https://en.wikipedia.org/wiki/SIMD
 *
 * There are different simd implementations (backends), which are auto-selected by seqan3::simd::simd_type_t.
 * Namely seqan3::detail::builtin_simd.
 *
 * \todo Make this public again. We made this documentation internal-documentation only for the 3.0.0 release, because
 * the API was not in shape yet. Remove /utility/simd/ exclusion from EXCLUDE_PATTERNS.
 *
 * \todo Describe more what SIMD is and how to use it.
 */

/*!\namespace seqan3::simd
 * \brief The SeqAn namespace for simd data types, algorithms and meta functions.
 * \sa https://en.wikipedia.org/wiki/SIMD What is SIMD conceptually?
 * \sa https://en.wikipedia.org/wiki/Streaming_SIMD_Extensions Which SIMD architectures exist?
 */

#pragma once

#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/simd/simd.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>
#include <seqan3/utility/simd/views/all.hpp>
