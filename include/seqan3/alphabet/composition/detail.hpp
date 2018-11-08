// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
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
// ==========================================================================

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides implementation detail for seqan3::union_composition and seqan3::cartesian_composition.
 */

#pragma once

#include <type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

// ------------------------------------------------------------------
// cartesian_composition_concept
// ------------------------------------------------------------------

/*!\interface seqan3::detail::cartesian_composition_concept <>
 * \extends seqan3::semi_alphabet_concept
 * \brief seqan3::cartesian_composition and its specialisations model this concept.
 * \ingroup alphabet
 *
 * \details
 *
 * This concept is necessary/helpful, because CRTP-specialisations cannot be tracked via regular inheritance or
 * specialisation mechanisms.
 */
//!\cond
template <typename t>
concept cartesian_composition_concept = requires
{
    typename t::seqan3_cartesian_components;
    typename t::seqan3_recursive_cartesian_components;
};
//!\endcond

// ------------------------------------------------------------------
// cartesian_components
// ------------------------------------------------------------------

/*!\brief Exposes for seqan3::cartesian_composition its components as a meta::list [base template].
 * \extends seqan3::TransformationTrait
 * \ingroup alphabet
 */
template <typename t>
struct cartesian_components;

/*!\brief Exposes for seqan3::cartesian_composition its components as a meta::list
 *        [specialisation for seqan3::cartesian_composition].
 * \extends seqan3::TransformationTrait
 * \ingroup alphabet
 */
template <cartesian_composition_concept t>
struct cartesian_components<t>
{
    //!\brief The returned type.
    using type = typename t::seqan3_cartesian_components;
};

// ------------------------------------------------------------------
// recursive_cartesian_components
// ------------------------------------------------------------------

/*!\brief Exposes for seqan3::cartesian_composition its components and those components' components (in the case of
 *        nested compositions) as a meta::list [base template].
 * \extends seqan3::TransformationTrait
 * \ingroup alphabet
 */
template <typename t>
struct recursive_cartesian_components;

/*!\brief Exposes for seqan3::cartesian_composition its components and those components' components (in the case of
 *        nested compositions) as a meta::list [specialisation for seqan3::cartesian_composition].
 * \extends seqan3::TransformationTrait
 * \ingroup alphabet
 */
template <cartesian_composition_concept t>
struct recursive_cartesian_components<t>
{
    //!\brief The returned type.
    using type = typename t::seqan3_recursive_cartesian_components;
};

// ------------------------------------------------------------------
// Callable concept helpers for meta::invoke
// ------------------------------------------------------------------

/*!\brief 'Callable' helper class that is invokable by meta::invoke.
 * Returns an std::true_type if the `type` is constructable from `T`.
 */
template <typename T>
struct constructible_from
{
    //!\brief The returned type when invoked.
    template <typename type>
    using invoke = std::integral_constant<bool, std::is_constructible_v<type, T>>;
};

/*!\brief 'Callable' helper class that is invokable by meta::invoke.
 * Returns an std::true_type if the `T` is implicitly convertible to `type`.
 */
template <typename T>
struct implicitly_convertible_from
{
    //!\brief The returned type when invoked.
    template <typename type>
    using invoke = std::integral_constant<bool, implicitly_convertible_to_concept<T, type>>;
};

/*!\brief 'Callable' helper class that is invokable by meta::invoke.
 * Returns an std::true_type if the `type` is assignable from `T`.
 */
template <typename T>
struct assignable_from
{
    //!\brief The returned type when invoked.
    template <typename type>
    using invoke = std::integral_constant<bool, weakly_assignable_concept<type, T>>;
};

/*!\brief 'Callable' helper class that is invokable by meta::invoke.
 * Returns an std::true_type if the `type` is weakly equality comparable to `T`.
 */
template <typename T>
struct weakly_equality_comparable_with
{
    //!\brief The returned type when invoked.
    template <typename type>
    using invoke = std::integral_constant<bool, std::detail::WeaklyEqualityComparableWith<type, T>>;
};

/*!\brief 'Callable' helper class that is invokable by meta::invoke.
 * Returns an std::true_type if the `type` is comparable via <,<=,>,>= to `T`.
 */
template <typename T>
struct weakly_ordered_with
{
    //!\brief The returned type when invoked.
    template <typename type>
    using invoke = std::integral_constant<bool, weakly_ordered_with_concept<type, T>>;
};

} // namespace seqan3::detail

// ------------------------------------------------------------------
// Forwards
// ------------------------------------------------------------------

namespace seqan3
{

// forward
template <typename ...alternative_types>
//!\cond
    requires (detail::constexpr_alphabet_concept<alternative_types> && ...) &&
             (sizeof...(alternative_types) >= 2)
             //TODO same char_type
//!\endcond
class union_composition;

template <typename derived_type,
          typename ...component_types>
//!\cond
    requires (detail::constexpr_semi_alphabet_concept<component_types> && ...)
//!\endcond
class cartesian_composition;

} // namespace seqan3
