// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <utility>

#include <seqan3/std/concepts>
#include <seqan3/std/type_traits>

namespace seqan3::detail
{

/*!\brief Applies a function to a range of arguments.
 * \ingroup algorithm
 *
 * \tparam unary_function_t The function type, like function pointers, functors and lambdas.
 * \tparam ...args_t        The parameter pack of the arguments. (Note that each argument type can be different)
 *
 * \param  fn               The function to call on every argument.
 * \param  args             The arguments as parameter pack.
 *
 * \details
 *
 * This function behaves like std::for_each but on parameter packs. The invocation(s) will be done without any loop.
 *
 * \include test/snippet/core/algorithm/for_each_value.cpp
 *
 * \attention The order of arguments is different to std::for_each. In std::for_each the first argument is a range and
 * the second one is a callable. Due to limitations with the ordering for parameter packs the order is the other way around.
 * the other way around.
 *
 * \sa https://en.cppreference.com/w/cpp/language/parameter_pack
 */
template <typename unary_function_t, typename ...args_t>
//!\cond
    requires (std::invocable<unary_function_t, args_t> && ... && true)
//!\endcond
constexpr void for_each_value(unary_function_t && fn, args_t && ...args)
{
    [[maybe_unused]] int r = (fn(std::forward<args_t>(args)), ..., 0);
}

/*!\brief Applies a function to a range of types.
 * \ingroup algorithm
 *
 * \tparam ...types         The parameter pack of types.
 * \tparam unary_function_t The function type, like function pointers, functors and lambdas.
 *
 * \param  fn               The function to call on every std::type_identity<type>.
 *
 * \details
 *
 * This function behaves like std::for_each but on types. The type will be wrapped into std::type_identity and passed as
 * argument. The invocation(s) will be done without any loop.
 *
 * This function can handle types which are incomplete, forward declared and not std::semiregular.
 *
 * \include test/snippet/core/algorithm/for_each_type.cpp
 *
 * \sa seqan3::detail::for_each_value
 */
template <typename ...types, typename unary_function_t>
//!\cond
    requires (std::invocable<unary_function_t, std::type_identity<types>> && ... && true)
//!\endcond
constexpr void for_each_type(unary_function_t && fn)
{
    for_each_value(std::forward<unary_function_t>(fn), std::type_identity<types>{}...);
}

/*!\brief Applies a function to a range of types contained in a list.
 * \ingroup algorithm
 *
 * \tparam type_list_t      The name of the list type.
 * \tparam ...types         The parameter pack of types within a list.
 * \tparam unary_function_t The function type, like function pointers, functors and lambdas.
 *
 * \param  fn               The function to call on every std::type_identity<type>.
 * \param  type_list        The list of types, i.e. seqan3::type_list, std::tuple. [unused]
 *
 * \details
 *
 * This function behaves like std::for_each but on types. The type will be wrapped into std::type_identity and passed as
 * argument. The invocation(s) will be done without any loop.
 *
 * This function can handle types which are incomplete, forward declared and not std::semiregular.
 *
 * \include test/snippet/core/algorithm/for_each_type_list.cpp
 *
 * \attention The order of arguments is different to std::for_each. In std::for_each the first argument is a range and
 * the second one is a callable. To make the interface consistent with seqan3::detail::for_each_value the order is changed.
 * order here, too.
 *
 * \sa seqan3::detail::for_each_value
 */
template <template <typename ...> typename type_list_t, typename ...args_t, typename unary_function_t>
//!\cond
    requires (std::invocable<unary_function_t, std::type_identity<args_t>> && ... && true)
//!\endcond
constexpr void for_each_type(unary_function_t && fn, type_list_t<args_t...> const & SEQAN3_DOXYGEN_ONLY(type_list))
{
    for_each_type<args_t...>(std::forward<unary_function_t>(fn));
}

} // namespace seqan3::detail
