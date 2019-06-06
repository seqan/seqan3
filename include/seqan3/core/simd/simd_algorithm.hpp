// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides algorithms to modify seqan3::simd::simd_type.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <utility>

#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_traits.hpp>

namespace seqan3::detail
{

//!\brief Helper function for seqan3::simd::fill.
//!\ingroup simd
template <Simd simd_t, size_t... I>
constexpr simd_t fill_impl(typename simd_traits<simd_t>::scalar_type const scalar, std::index_sequence<I...>)
{
    return simd_t{((void)I, scalar)...};
}

//!\brief Helper function for seqan3::simd::iota.
//!\ingroup simd
template <Simd simd_t, typename scalar_t, scalar_t... I>
constexpr simd_t iota_impl(scalar_t const offset, std::integer_sequence<scalar_t, I...>)
{
    return simd_t{static_cast<scalar_t>(offset + I)...};
}

} // namespace seqan3::detail

namespace seqan3
{

inline namespace simd
{

/*!\brief Fills a seqan3::simd::simd_type vector with a scalar value.
 * \tparam    simd_t The simd type which satisfies seqan3::simd::Simd.
 * \param[in] scalar The scalar value to fill the seqan3::simd::simd_type vector.
 * \ingroup simd
 *
 * \details
 *
 * \include test/snippet/core/simd/fill.cpp
 */
template <Simd simd_t>
constexpr simd_t fill(typename simd_traits<simd_t>::scalar_type const scalar)
{
    constexpr size_t length = simd_traits<simd_t>::length;
    return detail::fill_impl<simd_t>(scalar, std::make_index_sequence<length>{});
}

/*!\brief Fills a seqan3::simd::simd_type vector with the scalar values offset, offset+1, offset+2, ...
 * \tparam    simd_t The simd type which satisfies seqan3::simd::Simd.
 * \param[in] offset The scalar offset to fill the seqan3::simd::simd_type vector.
 * \ingroup simd
 *
 * \details
 *
 * \include test/snippet/core/simd/iota.cpp
 */
template <Simd simd_t>
constexpr simd_t iota(typename simd_traits<simd_t>::scalar_type const offset)
{
    constexpr size_t length = simd_traits<simd_t>::length;
    using scalar_type = typename simd_traits<simd_t>::scalar_type;
    return detail::iota_impl<simd_t>(offset, std::make_integer_sequence<scalar_type, length>{});
}

} // inline namespace simd

} // namespace seqan3
