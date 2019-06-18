// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Contains seqan3::simd::Simd
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{
//!\cond
// NOTE: this definition should be used for seqan3::Simd, but gcc has a bug that it will not fail silently if
// simd_t is a pointer to a incomplete type. Furthermore the is_pointer_v should prevent those cases by checking the type
// beforehand, but for some reasons the short-circuit semantic of `&&` does not work in this case and gcc still evaluates
// the requires clause which in turn triggers the error.
//
// If this concept is used directly on incomplete types it will produces this compiler error:
//     error: invalid use of incomplete type ‘struct incomplete::template_type<int>’
//          requires std::Same<decltype(a - b), simd_t>;
template <typename simd_t>
SEQAN3_CONCEPT Simd = requires (simd_t a, simd_t b)
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

} // namespace seqan3::detail

namespace seqan3
{

inline namespace simd
{

/*!\interface seqan3::simd::Simd <>
 * \brief The generic simd concept.
 * \ingroup simd
 *
 * \details
 *
 * seqan3::simd::Simd checks whether a given type is a simd type. One of the prerequisites is
 * that seqan3::simd::simd_traits is defined for this type.
 *
 * \if DEV
 * \todo Simplify concept to the seqan3::detail::Simd once gcc bug is fixed
 * \endif
 */
//!\cond
template <typename simd_t>
SEQAN3_CONCEPT Simd = !std::is_pointer_v<std::decay_t<simd_t>> && detail::Simd<simd_t>;
//!\endcond

} // inline namespace simd

} // namespace seqan3
