// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Meta-header for the simd module.
 */

#pragma once

/*!\defgroup simd Simd
 * \brief The simd module contains a unified interface to provide simd types and functions used in seqan3.
 * \see https://en.wikipedia.org/wiki/SIMD
 * \ingroup core
 *
 * There are different simd implementations (backends), which are auto-selected by seqan3::simd::simd_type_t.
 * \cond DEV
 * Namely seqan3::detail::builtin_simd.
 * \endcond
 *
 * \todo more details.
 */

#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/core/simd/simd_debug_stream.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/simd.hpp>

/*!\namespace seqan3::simd
 * \brief The SeqAn3 namespace for simd data types, algorithms and meta functions.
 * \sa https://en.wikipedia.org/wiki/SIMD What is SIMD conceptually?
 * \sa https://en.wikipedia.org/wiki/Streaming_SIMD_Extensions Which SIMD architectures exist?
 */
