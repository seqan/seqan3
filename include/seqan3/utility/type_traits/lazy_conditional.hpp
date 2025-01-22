// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides lazy template instantiation traits.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// lazy
// ----------------------------------------------------------------------------

/*!\brief An empty type whose only purpose is to hold an uninstantiated template plus its arguments.
 * \ingroup utility_type_traits
 * \tparam template_t The uninstantiated template.
 * \tparam spec_t     The arguments to template_t.
 */
template <template <typename...> typename template_t, typename... spec_t>
struct lazy
{};

// ----------------------------------------------------------------------------
// instantiate
// ----------------------------------------------------------------------------

/*!\brief A transformation trait that instantiates seqan3::lazy types. Base template is the identity transformation.
 * \ingroup utility_type_traits
 * \tparam t The type to operate on.
 * \implements seqan3::transformation_trait
 */
template <typename t>
struct instantiate : std::type_identity<t>
{};

/*!\brief A transformation trait that instantiates seqan3::lazy types.
 * \ingroup utility_type_traits
 * \tparam template_t The uninstantiated template.
 * \tparam spec_t     The arguments to template_t.
 * \implements seqan3::transformation_trait
 */
template <template <typename...> typename template_t, typename... spec_t>
struct instantiate<lazy<template_t, spec_t...>>
{
    //!\brief Return type of the trait [instantiates the template arguments].
    using type = template_t<spec_t...>;
};

/*!\brief A transformation trait that instantiates seqan3::lazy types. Transformation trait shortcut.
 * \tparam t The type to operate on.
 * \relates seqan3::detail::instantiate
 */
template <typename t>
    requires requires { typename instantiate<t>::type; }
using instantiate_t = typename instantiate<t>::type;

// ----------------------------------------------------------------------------
// instantiate_if
// ----------------------------------------------------------------------------

/*!\brief A transformation trait that instantiates seqan3::lazy types given a boolean condition.
 *        Base template is std::false_type.
 * \ingroup utility_type_traits
 * \tparam t The type to operate on.
 * \implements seqan3::transformation_trait
 */
template <typename t, bool condition>
struct instantiate_if : std::false_type
{};

/*!\brief A transformation trait that instantiates seqan3::lazy types given a boolean condition.
 *        If condition is true and parameter is not lazy, the type identity.
 * \ingroup utility_type_traits
 * \tparam t The type to operate on.
 * \implements seqan3::transformation_trait
 */
template <typename t>
struct instantiate_if<t, true> : std::type_identity<t>
{};

/*!\brief A transformation trait that instantiates seqan3::lazy types given a boolean condition.
 *        If condition is true and parameter is lazy, the instantiated type.
 * \ingroup utility_type_traits
 * \tparam template_t The uninstantiated template.
 * \tparam spec_t     The arguments to template_t.
 * \implements seqan3::transformation_trait
 */
template <template <typename...> typename template_t, typename... spec_t>
struct instantiate_if<lazy<template_t, spec_t...>, true>
{
    //!\brief Return type of the trait [instantiates the template arguments].
    using type = template_t<spec_t...>;
};

/*!\brief A transformation trait that instantiates seqan3::lazy types, conditionally. Transformation trait shortcut.
 * \tparam t The type to operate on.
 * \relates seqan3::detail::instantiate_if
 */
template <typename t, bool condition>
    requires requires { typename instantiate_if<t, condition>::type; }
using instantiate_if_t = typename instantiate_if<t, condition>::type;

/*!\brief A transformation trait that instantiates seqan3::lazy types, conditionally. Type trait shortcut.
 * \tparam t The type to operate on.
 * \relates seqan3::detail::instantiate_if
 */
template <typename t, bool condition>
    requires requires { instantiate_if_t<t, condition>::value; }
inline constexpr auto instantiate_if_v = instantiate_if_t<t, condition>::value;

// ----------------------------------------------------------------------------
// lazy_conditional
// ----------------------------------------------------------------------------

/*!\brief Behaves like std::conditional, but instantiates types wrapped in seqan3::lazy.
 * \ingroup utility_type_traits
 * \tparam decision   Whether to resolve to the first type or the second.
 * \tparam on_true_t  The return type in case `decision` is true.
 * \tparam on_false_t The return type in case `decision` is false.
 * \implements seqan3::transformation_trait
 *
 * \details
 *
 * ### Example
 *
 * \include test/snippet/utility/type_traits/lazy_conditional.cpp complete
 */
template <bool decision, typename on_true_t, typename on_false_t>
struct lazy_conditional : instantiate<std::conditional_t<decision, on_true_t, on_false_t>>
{};

/*!\brief Behaves like std::conditional_t, but instantiates types wrapped in seqan3::lazy. Transformation trait shortcut.
 * \tparam decision   Whether to resolve to the first type or the second.
 * \tparam on_true_t  The return type in case `decision` is true.
 * \tparam on_false_t The return type in case `decision` is false.
 * \relates seqan3::detail::lazy_conditional
 */
template <bool decision, typename on_true_t, typename on_false_t>
    requires requires { typename instantiate_t<std::conditional_t<decision, on_true_t, on_false_t>>; }
using lazy_conditional_t = instantiate_t<std::conditional_t<decision, on_true_t, on_false_t>>;

} // namespace seqan3::detail
