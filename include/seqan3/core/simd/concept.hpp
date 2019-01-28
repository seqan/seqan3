// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Contains seqan3::simd::simd_concept
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{

inline namespace simd
{

/*!\interface seqan3::simd::simd_concept <>
 * \brief The generic simd concept.
 * \ingroup simd
 *
 * \details
 *
 * seqan3::simd::simd_concept checks whether a given type is a simd type. One of the prerequisites is
 * that seqan3::simd::simd_traits is defined for this type.
 */
//!\cond
template <typename simd_t>
SEQAN3_CONCEPT simd_concept = requires (simd_t a, simd_t b)
{
    typename simd_traits<simd_t>::scalar_type;
    typename simd_traits<simd_t>::mask_type;
    typename simd_traits<simd_t>::swizzle_type;

    // require that static member variables are defined
    requires std::Integral<decltype(simd_traits<simd_t>::length)>;
    requires std::Integral<decltype(simd_traits<simd_t>::max_length)>;

    // assume array access that returns a scalar_type type
    { a[0] } -> typename simd_traits<simd_t>::scalar_type;

    // require comparison operators
    requires std::Same<decltype(a == b), typename simd_traits<simd_t>::mask_type>;
    requires std::Same<decltype(a != b), typename simd_traits<simd_t>::mask_type>;
    requires std::Same<decltype(a <  b), typename simd_traits<simd_t>::mask_type>;
    requires std::Same<decltype(a >  b), typename simd_traits<simd_t>::mask_type>;
    requires std::Same<decltype(a <= b), typename simd_traits<simd_t>::mask_type>;
    requires std::Same<decltype(a >= b), typename simd_traits<simd_t>::mask_type>;

    // require arithmetic operators
    requires std::Same<decltype(a + b), simd_t>;
    requires std::Same<decltype(a - b), simd_t>;
    requires std::Same<decltype(a * b), simd_t>;
    requires std::Same<decltype(a / b), simd_t>;
    requires std::Same<decltype(a += b), simd_t &>;
    requires std::Same<decltype(a -= b), simd_t &>;
    requires std::Same<decltype(a *= b), simd_t &>;
    requires std::Same<decltype(a /= b), simd_t &>;
};
//!\endcond

} // inline namespace simd

} // namespace seqan3
