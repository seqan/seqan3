// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::default_simd_length and seqan3::detail::default_simd_max_length
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/detail/int_types.hpp>

namespace seqan3::detail
{

/*!\brief seqan3 auto-detects the maximum number of packable `[u]int8_t` types.
 * \ingroup simd
 * \tparam simd_backend_t The name of the simd backend.
 *
 * \include test/snippet/core/simd/detail/default_simd_max_length.cpp
 *
 * This value is influenced by compiler flags like `-march=native`, `-msse4`,
 * `-mavx2`, etc and sets this value accordingly.
 * \sa https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html#x86-Options
 * for simd instruction sets and their flags.
 */
template <template <typename, size_t> typename simd_backend_t>
constexpr auto default_simd_max_length = 0u;

/*!\brief seqan3::detail::default_simd_length returns the default *length* depending
 * on the given *scalar_t* type, which is used in seqan3::simd::simd_type.
 * \ingroup simd
 * \tparam scalar_t The underlying type of a simd vector
 * \tparam simd_backend_t The name of the simd backend.
 */
template <typename scalar_t, template <typename, size_t> typename simd_backend_t>
constexpr auto default_simd_length = []
{
    constexpr auto max_length = default_simd_max_length<simd_backend_t>;
    if constexpr (max_length == 0)
        return min_viable_uint_v<1u>;
    else
        return min_viable_uint_v<max_length / sizeof(scalar_t)>;
}();

} // namespace seqan3::detail
