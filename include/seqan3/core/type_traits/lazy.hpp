// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides lazy template instantiation traits.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <seqan3/core/type_traits/function.hpp>
#include <seqan3/std/type_traits>

namespace seqan3::detail
{

/*!\addtogroup type_traits
 * \{
 */

// ----------------------------------------------------------------------------
// lazy
// ----------------------------------------------------------------------------

/*!\brief An empty type whose only purpose is to hold an uninstantiated template plus its arguments.
 * \tparam template_t The uninstantiated template.
 * \tparam spec_t     The arguments to template_t.
 */
template <template <typename ...> typename template_t, typename ...spec_t>
struct lazy
{};

// ----------------------------------------------------------------------------
// instantiate
// ----------------------------------------------------------------------------

/*!\brief A transformation trait that instantiates seqan3::lazy types. Base template is the identity transformation.
 * \tparam t The type to operate on.
 * \implements seqan3::transformation_trait
 */
template <typename t>
struct instantiate : std::type_identity<t>
{};

/*!\brief A transformation trait that instantiates seqan3::lazy types.
 * \tparam template_t The uninstantiated template.
 * \tparam spec_t     The arguments to template_t.
 * \implements seqan3::transformation_trait
 */
template <template <typename ...> typename template_t, typename ...spec_t>
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
//!\cond
    requires requires { typename instantiate<t>::type; }
//!\endcond
using instantiate_t = typename instantiate<t>::type;

// ----------------------------------------------------------------------------
// lazy_conditional
// ----------------------------------------------------------------------------

/*!\brief Behaves like std::conditional, but instantiates types wrapped in seqan3::lazy.
 * \tparam decision   Whether to resolve to the first type or the second.
 * \tparam on_true_t  The return type in case `decision` is true.
 * \tparam on_false_t The return type in case `decision` is false.
 * \implements seqan3::transformation_trait
 *
 * \details
 *
 * ### Example
 *
 * \include test/snippet/core/type_traits/lazy.cpp
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
//!\cond
    requires requires { typename instantiate_t<std::conditional_t<decision, on_true_t, on_false_t>>; }
//!\endcond
using lazy_conditional_t = instantiate_t<std::conditional_t<decision, on_true_t, on_false_t>>;

/*!\brief An unary type trait that tests whether a template class can be instantiated with the given template type
 *        parameters.
 * \implements seqan3::unary_type_trait
 * \tparam query_t The type of the template class to test.
 * \tparam args_t  The template parameter pack to instantiate the template class with.
 *
 * \details
 *
 * Note, this unary type trait can be used in a seqan3::detail::lazy_conditional expression to check if instantiating
 * a template class with specific template arguments would result in a valid template definition. Thus, the template
 * parameters of the checked class must be constrained accordingly.
 *
 * ### Example
 *
 * \include test/snippet/core/type_traits/is_instantiable_with.cpp
 */
template <template <typename ...> typename query_t, typename ...args_t>
struct is_instantiable_with :
//!\cond
    public std::false_type
//!\endcond
{};

//!\cond
template <template <typename ...> typename query_t, typename ...args_t>
    requires requires { typename std::type_identity<query_t<args_t...>>::type; }
struct is_instantiable_with<query_t, args_t...> : public std::true_type
{};
//!\endcond

/*!\brief Helper variable template for seqan3::detail::is_instantiable_with.
 * \tparam query_t The type of the template class to test.
 * \tparam args_t  The template parameter pack to instantiate the template class with.
 * \relates seqan3::detail::is_instantiable_with
 */
template <template <typename ...> typename query_t, typename ...args_t>
inline constexpr bool is_instantiable_with_v = is_instantiable_with<query_t, args_t...>::value;
//!\}
} // namespace seqan3::detail
