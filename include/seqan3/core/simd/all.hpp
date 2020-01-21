// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Meta-header for the simd module.
 */

#pragma once

//!\cond DEV

/*!\defgroup simd Simd
 * \brief The simd module contains a unified interface to provide simd types and functions used in seqan3.
 * \see https://en.wikipedia.org/wiki/SIMD
 * \ingroup core
 *
 * \todo Make this public again. We made this documentation internal-documentation only for the 3.0.0 release, because
 * the API was not in shape yet. Remove the `cond DEV` and `todo` commands and remove `seqan3::simd` from
 * SEQAN3_DOXYGEN_EXCLUDE_SYMBOLS in `seqan3-doxygen.cmake`.
 *
 * There are different simd implementations (backends), which are auto-selected by seqan3::simd::simd_type_t.
 * \if DEV
 * Namely seqan3::detail::builtin_simd.
 * \endif
 *
 * \todo Describe more what SIMD is and how to use it.
 */

#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/debug_stream_simd.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/view_to_simd.hpp>

/*!\namespace seqan3::simd
 * \brief The SeqAn namespace for simd data types, algorithms and meta functions.
 * \sa https://en.wikipedia.org/wiki/SIMD What is SIMD conceptually?
 * \sa https://en.wikipedia.org/wiki/stream_REMOVEMEing_SIMD_Extensions Which SIMD architectures exist?
 */
//!\endcond
