// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides implementation detail for seqan3::alphabet_variant and seqan3::alphabet_tuple_base.
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/type_list/traits.hpp>
#include <seqan3/core/type_traits/lazy.hpp>

namespace seqan3::detail
{

// ------------------------------------------------------------------
// alphabet_tuple_base_specialisation
// ------------------------------------------------------------------

/*!\interface seqan3::detail::alphabet_tuple_base_specialisation <>
 * \brief seqan3::alphabet_tuple_base and its derivates model this concept.
 * \ingroup composite
 *
 * \details
 *
 * This concept is necessary/helpful, because CRTP-specialisations cannot easily be tracked via regular inheritance or
 * specialisation mechanisms.
 */
//!\cond
template <typename t>
SEQAN3_CONCEPT alphabet_tuple_base_specialisation = requires
{
    requires t::seqan3_alphabet_tuple_base_specialisation;
};
//!\endcond

// ------------------------------------------------------------------
// required_types
// ------------------------------------------------------------------

/*!\brief A seqan3::type_list with types that the given type depends on.
 * \implements seqan3::transformation_trait
 * \ingroup composite
 *
 * \details
 *
 * The list is empty by default. This trait maybe used in metaprogramming to indicate that certain types need to be
 * complete and not depend on the given type to avoid recursive template instantiation.
 */
template <typename t>
struct required_types
{
    //!\brief The returned type.
    using type = type_list<>;
};

/*!\brief A seqan3::type_list with types that the given type depends on.
 *        [specialisation for seqan3::alphabet_variant and derivates of seqan3::alphabet_tuple_base].
 * \implements seqan3::transformation_trait
 * \ingroup composite
 *
 * \details
 *
 * Exposes for seqan3::alphabet_tuple_base its components and for seqan3::alphabet_variant its alternatives.
 */
template <typename t>
//!\cond
    requires requires { typename t::seqan3_required_types; }
//!\endcond
struct required_types<t>
{
    //!\brief The returned type.
    using type = typename t::seqan3_required_types;
};

/*!\brief A seqan3::type_list with types that the given type depends on. [Trait shortcut]
 * \relates seqan3::detail::required_types
 */
template <typename t>
using required_types_t = typename required_types<t>::type;

// ------------------------------------------------------------------
// recursive_required_types
// ------------------------------------------------------------------

//TODO: This can be replaced with metaprogramming magic once a few more functions land in list_traits.

/*!\brief Like seqan3::detail::required_types, but recursive.
 * \implements seqan3::transformation_trait
 * \ingroup composite
 */
template <typename t>
struct recursive_required_types
{
    //!\brief The returned type.
    using type = type_list<>;
};

/*!\brief Like seqan3::detail::required_types, but recursive.
 * \implements seqan3::transformation_trait
 * \ingroup composite
 */
template <typename t>
//!\cond
    requires requires
    {
        typename t::seqan3_recursive_required_types;
    }
//!\endcond
struct recursive_required_types<t>
{
    //!\brief The returned type.
    using type = typename t::seqan3_recursive_required_types;
};

/*!\brief Shortcut for seqan3::detail::recursive_required_types.
 * \relates seqan3::detail::recursive_required_types
 */
template <typename t>
using recursive_required_types_t = typename recursive_required_types<t>::type;

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
    using invoke = std::integral_constant<bool, implicitly_convertible_to<T, type>>;
};

/*!\brief 'Callable' helper class that is invokable by meta::invoke.
 * Returns an std::true_type if the `type` is assignable from `T`.
 */
template <typename T>
struct assignable_from
{
    //!\brief The returned type when invoked.
    template <typename type>
    using invoke = std::integral_constant<bool, weakly_assignable_from<type, T>>;
};

/*!\brief 'Callable' helper class that is invokable by meta::invoke.
 * Returns an std::true_type if the `type` is weakly equality comparable to `T`.
 */
template <typename T>
struct weakly_equality_comparable_with_
{
    //!\brief The returned type when invoked.
    template <typename type>
    using invoke = std::integral_constant<bool, weakly_equality_comparable_with<type, T>>;
};

/*!\brief 'Callable' helper class that is invokable by meta::invoke.
 * Returns an std::true_type if the `type` is comparable via <,<=,>,>= to `T`.
 */
template <typename T>
struct weakly_ordered_with_
{
    //!\brief The returned type when invoked.
    template <typename type>
    using invoke = std::integral_constant<bool, weakly_ordered_with<type, T>>;
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
    requires (detail::writable_constexpr_alphabet<alternative_types> && ...) &&
             (std::regular<alternative_types> && ...) &&
             (sizeof...(alternative_types) >= 2)
             //TODO same char_type
//!\endcond
class alphabet_variant;

template <typename derived_type,
          typename ...component_types>
//!\cond
    requires (detail::writable_constexpr_semialphabet<component_types> && ...) &&
             (std::regular<component_types> && ...)
//!\endcond
class alphabet_tuple_base;

} // namespace seqan3
