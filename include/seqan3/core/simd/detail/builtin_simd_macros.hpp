// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

//!\cond DEV

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides #SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION
 */

#pragma once

#include <seqan3/core/platform.hpp>

/*!\brief #SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_CALL calls the
 * given macro *class_definition* with parameters needed to define
 * [vector extension](https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html)
 * types.
 * \ingroup simd
 * \param  class_definition the macro to call
 * \param  scalar_t   same as seqan3::simd_traits<builtin_simd_t>::scalar_t
 * \param  max_length same as seqan3::simd_traits<builtin_simd_t>::max_length
 *
 * Calls *class_definition* with the following parameters:
 * \li *scalar_t*   same as seqan3::simd_traits<builtin_simd_t>::scalar_t
 * \li *length*     same as seqan3::simd_traits<builtin_simd_t>::length
 * \li *max_length* same as seqan3::simd_traits<builtin_simd_t>::max_length
 * \li *simd_t*     the type of the [vector extension]
 *                  (https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html),
 *                  e.g. `scalar_t __attribute__ ((__vector_size__(max_length)))`
 *
 * \par Called by
 * #SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_BY_MAX_LENGTH
 */
#define SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_CALL(class_definition, scalar_t, max_length)/*
*/class_definition(scalar_t, max_length/sizeof(scalar_t), max_length, scalar_t __attribute__ ((__vector_size__(max_length))))

/*!\brief #SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_BY_MAX_LENGTH
 * will call #SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_CALL for all
 * `[u]intX_t` types for a given *max_length*
 * \ingroup simd
 * \param  class_definition the macro to call
 * \param  max_length same as seqan3::simd_traits<builtin_simd_t>::max_length
 *
 * \par Called by
 * #SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION
 */
#define SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_BY_MAX_LENGTH(class_definition, max_length) \
    SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_CALL(class_definition, int8_t  , max_length) \
    SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_CALL(class_definition, uint8_t , max_length) \
    SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_CALL(class_definition, int16_t , max_length) \
    SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_CALL(class_definition, uint16_t, max_length) \
    SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_CALL(class_definition, int32_t , max_length) \
    SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_CALL(class_definition, uint32_t, max_length) \
    SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_CALL(class_definition, int64_t , max_length) \
    SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_CALL(class_definition, uint64_t, max_length)

/*!\brief #SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION
 * will call *class_definition* for all `[u]intX_t` types and supported simd
 * extensions (e.g. sse4, avx2, avx512).
 * \ingroup simd
 * \param  class_definition the macro to call
 *
 * Calls *class_definition* with the following parameters:
 * \li *scalar_t*   same as seqan3::simd_traits<builtin_simd_t>::scalar_t
 * \li *length*     same as seqan3::simd_traits<builtin_simd_t>::length
 * \li *max_*length same as seqan3::simd_traits<builtin_simd_t>::max_length
 * \li *simd_t*     the type of the [vector extension]
 * (https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html),
 * e.g. `scalar_t __attribute__ ((__vector_size__(max_length)))`
 *
 * \sa seqan3::detail::builtin_simd for an example
 */
#define SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION(class_definition) \
    SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_BY_MAX_LENGTH(class_definition, 16) /*128bit = 16 * 8bit*/ \
    SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_BY_MAX_LENGTH(class_definition, 32) /*256bit = 32 * 8bit*/ \
    SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION_BY_MAX_LENGTH(class_definition, 64) /*512bit = 64 * 8bit*/

//!\endcond
