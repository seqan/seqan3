// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides algorithms for meta programming, parameter packs and seqan3::type_list.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <type_traits>
#include <utility>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

//-----------------------------------------------------------------------------
// all_of
//-----------------------------------------------------------------------------

/*!\brief Tests whether a given predicate evaluates to `true` for each element in the function parameter pack.
 * \ingroup utility_type_pack
 *
 * \tparam unary_predicate_t The function type, like function pointers, functors and lambdas;
 *                           must model std::predicate expanded on each argument type.
 * \tparam pack_t The parameter pack of the arguments (each argument type can be different).
 *
 * \param[in] fn The predicate to evaluate for every argument.
 * \param[in] args The parameter pack.
 *
 * \returns `true` if the predicate returns `true` for each type in the type list, `false` otherwise.
 *
 * \details
 *
 * This function behaves like std::all_of but on parameter packs. The invocation(s) will be done without any loop.
 *
 * ### Example
 *
 * \include test/snippet/utility/type_pack/detail/type_pack_algorithm_all_of.cpp
 *
 * ### Complexity
 *
 * Linear in the number of elements in the pack.
 *
 * \attention Opposed to the std::all_of the argument order is changed, such that the first argument is the unary
 * predicate to invoke on each argument followed by the arguments.
 * This is due to a constraint in the c++ language regarding parameter packs.
 *
 * \sa https://en.cppreference.com/w/cpp/language/parameter_pack
 */
template <typename unary_predicate_t, typename... pack_t>
    requires (std::predicate<unary_predicate_t, pack_t> && ...)
constexpr bool all_of(unary_predicate_t && fn, pack_t &&... args)
{
    return (fn(std::forward<pack_t>(args)) && ...);
}

//-----------------------------------------------------------------------------
// for_each
//-----------------------------------------------------------------------------

/*!\brief Applies a function to each element of the given function parameter pack.
 * \ingroup utility_type_pack
 *
 * \tparam unary_function_t The function type, like function pointers, functors and lambdas.
 * \tparam pack_t The parameter pack of the arguments (each argument type can be different).
 *
 * \param[in] fn The function to call on every argument.
 * \param[in] args The parameter pack.
 *
 * \details
 *
 * This function behaves like std::for_each but on parameter packs. The invocation(s) will be done without any loop.
 *
 * ### Example
 *
 * \include test/snippet/utility/type_pack/detail/type_pack_algorithm_for_each.cpp
 *
 * ### Complexity
 *
 * Linear in the number of elements in the pack.
 *
 * \attention Opposed to the std::for_each the argument order is changed, such that the first argument is the unary
 * function to invoke on each argument followed by the arguments.
 * This is due to a constraint in the c++ language regarding parameter packs.
 *
 * \sa https://en.cppreference.com/w/cpp/language/parameter_pack
 */
template <typename unary_function_t, typename... pack_t>
    requires (std::invocable<unary_function_t, pack_t> && ...)
constexpr void for_each(unary_function_t && fn, pack_t &&... args)
{
    (fn(std::forward<pack_t>(args)), ...);
}

} // namespace seqan3::detail
