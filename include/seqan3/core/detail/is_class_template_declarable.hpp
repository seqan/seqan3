// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a type trait for verifying valid template declarations.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\addtogroup core
 * \{
 */

/*!\brief An unary type trait that tests whether a template class can be declared with the given template type
 *        parameters.
 * \implements seqan3::unary_type_trait
 * \tparam query_t The type of the template class to test.
 * \tparam args_t  The template parameter pack to instantiate the template class with.
 *
 * \details
 *
 * Note, this unary type trait can be used in a seqan3::detail::lazy_conditional expression to check if instantiating
 * a template class with specific template arguments would result in a valid template declaration. Thus, the template
 * parameters of the checked class must be constrained accordingly. It is, however, not possible to test if the
 * the resulting type is incomplete or not, such that it can not be tested if an instance of the class template with
 * the given template arguments can be actually created.
 *
 * ### Example
 *
 * \include test/snippet/core/detail/is_class_template_declarable_with.cpp
 */
template <template <typename ...> typename query_t, typename ...args_t>
struct is_class_template_declarable_with :
//!\cond
    public std::false_type
//!\endcond
{};

//!\cond
template <template <typename ...> typename query_t, typename ...args_t>
requires requires { typename std::type_identity<query_t<args_t...>>::type; }
struct is_class_template_declarable_with<query_t, args_t...> : public std::true_type
{};
//!\endcond

/*!\brief Helper variable template for seqan3::detail::is_class_template_declarable_with.
 * \tparam query_t The type of the template class to test.
 * \tparam args_t  The template parameter pack to instantiate the template class with.
 * \relates seqan3::detail::is_class_template_declarable_with
 */
template <template <typename ...> typename query_t, typename ...args_t>
inline constexpr bool is_class_template_declarable_with_v =
    is_class_template_declarable_with<query_t, args_t...>::value;

//!\}

} // namespace seqan3::detail
