// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides algorithms for meta programming, parameter packs and seqan3::type_list.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <utility>

#include <seqan3/core/type_list/type_list.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/type_traits>

namespace seqan3::detail
{

//!\cond
template <typename type_list_t>
struct type_list_expander;
//!\endcond

/*!\brief Helper class to invoke a meta algorithm on the types contained in a seqan3::type_list.
 * \tparam type_list_t The type list given as a template template parameter.
 * \tparam args_t The template arguments contained in `type_list_t` to apply the target function on.
 *
 * \details
 *
 * The meta algorithms provide a parameter pack version and a seqan3::type_list version.
 * In the type list version the algorithm is called on the types contained in the enclosing type list after they have
 * been wrapped in std::type_identity. Thus, the type list version uses the parameter pack version to call the passed
 * function on the types.
 *
 * Concrete this helper does the following:
 *  * expand the types within the type list
 *  * wrap each type in std::identity and instantiate it (i.e. transform a type into a value)
 *  * call a target function with the instances of std::identity (as a pack).
 *
 * This is a technical trick to make a type representable as a value. Instantiating a type might not always work
 * because not every type provides a default constructor. In addition it is possible to use incomplete types as well.
 */
template <template <typename ...> typename type_list_t, typename ...args_t>
struct type_list_expander<type_list_t<args_t...>>
{
    /*!\brief Invokes the actual function by passing the types as instances of std::type_identity to the target
     *        function.
     * \tparam fn_t The type of the function to be called with the expanded list of types.
     * \param[in] fn The function to be called.
     * \returns The std::invoke_result of calling `fn` with the expanded list of types.
     *
     * \details
     *
     * Invokes `fn` by passing the expanded types wrapped in std::type_identity.
     */
    template <typename fn_t>
    //!\cond
        requires std::invocable<fn_t, std::type_identity<args_t>...>
    //!\endcond
    static constexpr std::invoke_result_t<fn_t, std::type_identity<args_t>...> invoke_on_type_identities(fn_t && fn)
    {
        return fn(std::type_identity<args_t>{}...);
    }
};

//-----------------------------------------------------------------------------
// all_of
//-----------------------------------------------------------------------------

/*!\brief Tests whether a given predicate evaluates to `true` for each element in the function parameter pack.
 * \ingroup core
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
 * \include test/snippet/core/detail/all_of_in_pack.cpp
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
template <typename unary_predicate_t, typename ...pack_t>
//!\cond
    requires (std::predicate<unary_predicate_t, pack_t> && ...)
//!\endcond
constexpr bool all_of(unary_predicate_t && fn, pack_t && ...args)
{
    return (fn(std::forward<pack_t>(args)) && ...);
}

/*!\brief Tests whether a given predicate evaluates to `true` for each type in a seqan3::type_list.
 * \ingroup core
 *
 * \tparam list_t A type list; must model seqan3::detail::type_list_specialisation
 * \tparam unary_predicate_t The function type, like function pointers, functors and lambdas;
 *                           must model std::predicate expanded on each argument type wrapped in std::type_identity.
 *
 * \param[in] fn The predicate called for every type in the seqan3::type_list.
 *
 * \returns `true` if the predicate returns `true` for each type in the type list, `false` otherwise.
 *
 * \details
 *
 * This function operates on types instead of values.
 * The following steps are performed to call the passed predicate on the types contained in the type list:
 *
 *  * expand the types within the type list
 *  * wrap each type in std::identity and instantiate it (i.e. transform a type into a value)
 *  * call the parameter pack version of seqan3::detail::all_of with the instances of std::identity (as a pack).
 *
 * Note that wrapping the types in std::type_identity is a technical trick to make a type representable as a value.
 * Instantiating a type might not work because they might not be std::default_constructible.
 * In addition it is possible, to invoke the predicate on incomplete types.
 *
 * ### Example
 *
 * \include test/snippet/core/detail/all_of_in_type_list.cpp
 *
 * ### Complexity
 *
 * Linear in the number of types in the seqan3::type_list.
 *
 * [Compile-time complexity: Linear number of template instantiations.]
 */
template <typename type_list_t, typename unary_predicate_t>
[[nodiscard]] constexpr bool all_of(unary_predicate_t && fn)
//!\cond
    requires type_list_specialisation<type_list_t>
//!\endcond
{
    return type_list_expander<type_list_t>::invoke_on_type_identities([&] (auto && ...type_identities)
    {
        return all_of(fn, std::forward<decltype(type_identities)>(type_identities)...);
    });
}

//-----------------------------------------------------------------------------
// for_each
//-----------------------------------------------------------------------------

/*!\brief Applies a function to each element of the given function parameter pack.
 * \ingroup core
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
 * \include test/snippet/core/detail/for_each_in_pack.cpp
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
template <typename unary_function_t, typename ...pack_t>
//!\cond
    requires (std::invocable<unary_function_t, pack_t> && ...)
//!\endcond
constexpr void for_each(unary_function_t && fn, pack_t && ...args)
{
    (fn(std::forward<pack_t>(args)), ...);
}

/*!\brief Applies a function element wise to all types of a type list.
 * \ingroup core
 *
 * \tparam list_t A type list; must model seqan3::detail::type_list_specialisation.
 * \tparam unary_function_t The function type, like function pointers, functors and lambdas; must model
 *                          std::invocable on each type of the type list wrapped in std::type_identity.
 *
 * \param[in] fn The function to call on every type contained in the list.
 *
 * \details
 *
 * This function operates on types instead of values.
 * The following steps are performed to call the passed unary function on the types contained in the type list:
 *
 *  * expand the types within the type list
 *  * wrap each type in std::identity and instantiate it (i.e. transform a type into a value)
 *  * call the parameter pack version of seqan3::detail::for_each with the instances of std::identity (as a pack).
 *
 * Note that wrapping the types in std::type_identity is a technical trick to make a type representable as a value.
 * Instantiating a type might not work because they might not be std::default_constructible.
 * In addition, it is possible to invoke the unary function on incomplete types.
 *
 * ### Example
 *
 * \include test/snippet/core/detail/for_each_in_type_list.cpp
 *
 * ### Complexity
 *
 * Linear in the number of types in the seqan3::type_list.
 *
 * [Compile-time complexity: Linear number of template instantiations.]
 *
 * \sa seqan3::detail::for_each
 */
template <typename type_list_t, typename unary_function_t>
//!\cond
    requires type_list_specialisation<type_list_t>
//!\endcond
constexpr void for_each(unary_function_t && fn)
{
    type_list_expander<type_list_t>::invoke_on_type_identities([&] (auto && ...type_identities)
    {
        for_each(fn, std::forward<decltype(type_identities)>(type_identities)...);
    });
}

}  // namespace seqan3::detail
