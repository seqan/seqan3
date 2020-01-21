// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::default_simd_backend
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/simd/detail/builtin_simd.hpp>

namespace seqan3::detail
{

/*!\brief seqan3::detail::default_simd_backend is the default used implementation of
 * seqan3::simd::simd_type.
 * \ingroup simd
 * \tparam scalar_t The underlying type of a simd vector
 * \tparam length The number of packed values in a simd vector
 */
template <typename scalar_t, size_t length>
using default_simd_backend = builtin_simd<scalar_t, length>;
} // namespace seqan3::detail
