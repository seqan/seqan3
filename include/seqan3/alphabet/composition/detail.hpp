// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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
SEQAN3_CONCEPT cartesian_composition_concept = requires
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
