// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::default_simd_backend
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/utility/simd/detail/builtin_simd.hpp>

namespace seqan3::detail
{

/*!\brief seqan3::detail::default_simd_backend is the default used implementation of
 * seqan3::simd::simd_type.
 * \ingroup utility_simd
 * \tparam scalar_t The underlying type of a simd vector
 * \tparam length The number of packed values in a simd vector
 */
template <typename scalar_t, size_t length>
using default_simd_backend = builtin_simd<scalar_t, length>;
} // namespace seqan3::detail
