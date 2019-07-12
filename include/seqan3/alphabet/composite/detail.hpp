// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides implementation detail for seqan3::alphabet_variant and seqan3::alphabet_tuple_base.
 */

#pragma once

#include <type_traits>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/core/concept/core_language.hpp>

namespace seqan3::detail
{

// ------------------------------------------------------------------
// AlphabetTupleBase
// ------------------------------------------------------------------

/*!\interface seqan3::detail::AlphabetTupleBase <>
 * \extends seqan3::Semialphabet
 * \brief seqan3::alphabet_tuple_base and its specialisations model this concept.
 * \ingroup alphabet
 *
 * \details
 *
 * This concept is necessary/helpful, because CRTP-specialisations cannot be tracked via regular inheritance or
 * specialisation mechanisms.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT AlphabetTupleBase = requires
{
    typename t::seqan3_tuple_components;
    typename t::seqan3_recursive_tuple_components;
};
//!\endcond

// ------------------------------------------------------------------
// tuple_components
// ------------------------------------------------------------------

/*!\brief Exposes for seqan3::alphabet_tuple_base its components as a meta::list.
 * \implements seqan3::TransformationTrait
 * \ingroup alphabet
 */
template <typename t>
struct tuple_components;

/*!\brief Exposes for seqan3::alphabet_tuple_base its components as a meta::list
 *        [specialisation for seqan3::alphabet_tuple_base].
 * \implements seqan3::TransformationTrait
 * \ingroup alphabet
 */
template <AlphabetTupleBase t>
struct tuple_components<t>
{
    //!\brief The returned type.
    using type = typename t::seqan3_tuple_components;
};

// ------------------------------------------------------------------
// recursive_tuple_components
// ------------------------------------------------------------------

/*!\brief Exposes for seqan3::alphabet_tuple_base its components and those components' components (in the case of
 *        nested composites) as a meta::list [base template].
 * \implements seqan3::TransformationTrait
 * \ingroup alphabet
 */
template <typename t>
struct recursive_tuple_components;

/*!\brief Exposes for seqan3::alphabet_tuple_base its components and those components' components (in the case of
 *        nested composites) as a meta::list [specialisation for seqan3::alphabet_tuple_base].
 * \implements seqan3::TransformationTrait
 * \ingroup alphabet
 */
template <AlphabetTupleBase t>
struct recursive_tuple_components<t>
{
    //!\brief The returned type.
    using type = typename t::seqan3_recursive_tuple_components;
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
    using invoke = std::integral_constant<bool, ImplicitlyConvertibleTo<T, type>>;
};

/*!\brief 'Callable' helper class that is invokable by meta::invoke.
 * Returns an std::true_type if the `type` is assignable from `T`.
 */
template <typename T>
struct assignable_from
{
    //!\brief The returned type when invoked.
    template <typename type>
    using invoke = std::integral_constant<bool, WeaklyAssignable<type, T>>;
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
    using invoke = std::integral_constant<bool, WeaklyOrderedWith<type, T>>;
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
    requires (detail::WritableConstexprAlphabet<alternative_types> && ...) &&
             (!std::is_reference_v<alternative_types> && ...) &&
             (sizeof...(alternative_types) >= 2)
             //TODO same char_type
//!\endcond
class alphabet_variant;

template <typename derived_type,
          typename ...component_types>
//!\cond
    requires (detail::WritableConstexprSemialphabet<component_types> && ...) &&
             (!std::is_reference_v<component_types> && ...)
//!\endcond
class alphabet_tuple_base;

} // namespace seqan3
