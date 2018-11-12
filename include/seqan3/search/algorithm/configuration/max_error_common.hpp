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

/*!\file
 * \brief Provides the error types for maximum number of errors.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/std/concepts>

namespace seqan3::search_cfg
{

/*!\brief A strong type of underlying type `uint8_t` or `double` that represents the number or rate of total errors.
 * \tparam value_t The underlying type
 * \ingroup search_configuration
 */
template <typename value_t>
//!\cond
    requires Arithmetic<value_t>
//!\endcond
struct total : detail::strong_type<value_t, total<value_t>, detail::strong_type_skill::convert>
{
    using detail::strong_type<value_t, total<value_t>, detail::strong_type_skill::convert>::strong_type;
};

/*!\name Template argument type deduction guides
 * \relates seqan3::search_cfg::total
 * \{
 */
template <std::Integral value_t>
total(value_t) -> total<uint8_t>;

template <typename value_t>
//!\cond
    requires FloatingPoint<value_t>
//!\endcond
total(value_t) -> total<double>;
//!\}

/*!\brief A strong type of underlying type `uint8_t` or `double` that represents the number or rate of
 *        substitutions.
 * \tparam value_t The underlying type
 * \ingroup search_configuration
 */
template <typename value_t>
//!\cond
    requires Arithmetic<value_t>
//!\endcond
struct substitution : detail::strong_type<value_t, substitution<value_t>, detail::strong_type_skill::convert>
{
    using detail::strong_type<value_t, substitution<value_t>, detail::strong_type_skill::convert>::strong_type;
};

/*!\name Template argument type deduction guides
 * \relates seqan3::search_cfg::substitution
 * \{
 */
template <std::Integral value_t>
substitution(value_t) -> substitution<uint8_t>;

template <typename value_t>
//!\cond
    requires FloatingPoint<value_t>
//!\endcond
substitution(value_t) -> substitution<double>;
//!\}

/*!\brief A strong type of underlying type `uint8_t` or `double` that represents the number or rate of insertions.
 * \tparam value_t The underlying type
 * \ingroup search_configuration
 */
template <typename value_t>
//!\cond
    requires Arithmetic<value_t>
//!\endcond
struct insertion : detail::strong_type<value_t, insertion<value_t>, detail::strong_type_skill::convert>
{
    using detail::strong_type<value_t, insertion<value_t>, detail::strong_type_skill::convert>::strong_type;
};

/*!\name Template argument type deduction guides
 * \relates seqan3::search_cfg::total
 * \{
 */
template <std::Integral value_t>
insertion(value_t) -> insertion<uint8_t>;

template <typename value_t>
//!\cond
    requires FloatingPoint<value_t>
//!\endcond
insertion(value_t) -> insertion<double>;
//!\}

/*!\brief A strong type of underlying type `uint8_t` or `double` that represents the number or rate of deletions.
 * \tparam value_t The underlying type
 * \ingroup search_configuration
 */
template <typename value_t>
//!\cond
    requires Arithmetic<value_t>
//!\endcond
struct deletion : detail::strong_type<value_t, deletion<value_t>, detail::strong_type_skill::convert>
{
    using detail::strong_type<value_t, deletion<value_t>, detail::strong_type_skill::convert>::strong_type;
};

/*!\name Template argument type deduction guides
 * \relates seqan3::search_cfg::deletion
 * \{
 */
template <std::Integral value_t>
deletion(value_t) -> deletion<uint8_t>;

template <typename value_t>
//!\cond
    requires FloatingPoint<value_t>
//!\endcond
deletion(value_t) -> deletion<double>;
//!\}

} // namespace seqan3::search_cfg
