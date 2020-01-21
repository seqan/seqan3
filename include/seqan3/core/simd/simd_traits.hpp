// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::simd::simd_traits
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{

inline namespace simd
{

/*!\brief seqan3::simd::simd_traits is the trait class that provides uniform interface
 * to the properties of simd_t types.
 * \ingroup simd
 * \tparam simd_t The simd type that satisfies seqan3::simd::simd_concept.
 *
 * The class defines the following member variables and types:
 * * scalar_type - the underlying type of a simd vector
 * * length - the number of packed values in a simd vector
 * * max_length - the maximum number of packable values in a simd vector, if the underlying type would be [u]int8_t
 * * mask_type - the type returned by comparison operators
 * * swizzle_type - the type used to define how to swizzle a simd vector
 *
 * \include test/snippet/core/simd/simd_traits.cpp
 */
template <typename simd_t>
struct simd_traits
#if SEQAN3_DOXYGEN_ONLY(1)0
{
    /*!\brief The underlying type of a simd vector (is not defined if *simd_t*
     * does not model *seqan3::simd::simd*)
     */
    using scalar_type = IMPLEMENTATION_DEFINED;
    /*!\brief The number of packed values in a simd vector (is not defined if
     * *simd_t* does not model *seqan3::simd::simd*)
     */
    static constexpr auto length = IMPLEMENTATION_DEFINED;
    /*!\brief The maximum number of packable values in a simd vector, if the
     * underlying type would be *[u]int8_t* (is not defined if *simd_t* does not
     * model *seqan3::simd::simd*)
     */
    static constexpr auto max_length = IMPLEMENTATION_DEFINED;
    /*!\brief The type returned by comparison operators (is not defined if
     * *simd_t* does not model *seqan3::simd::simd*)
     */
    using mask_type = IMPLEMENTATION_DEFINED;
    /*!\brief The type used to define how to swizzle a simd vector (is not
     * defined if *simd_t* does not model *seqan3::simd::simd*)
     */
    using swizzle_type = IMPLEMENTATION_DEFINED;
}
#endif
;

} // inline namespace simd

} // namespace seqan3
