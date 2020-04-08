// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various type traits for use on functions.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <functional>

#include <seqan3/core/type_list/traits.hpp>

namespace seqan3
{
// ----------------------------------------------------------------------------
// function_traits
// ----------------------------------------------------------------------------

//!\cond
template <typename function_t>
struct function_traits;
//!\endcond

/*!\brief A traits class to provide a uniform interface to the properties of a function type.
 * \ingroup type_traits
 * \tparam return_t The return type of the function.
 * \tparam args_t A template parameter pack over the argument types of the function.
 *
 * \details
 *
 * seqan3::function_traits is the type trait class that provides a uniform interface to the properties of
 * a std::function type, a lambda type, a function type or a function pointer type.
 * This makes it possible to access the return type and the argument types of the stored target function.
 * The function types must be complete, i.e. all argument types and the return type must be known, otherwise
 * this traits class is incomplete.
 *
 * ### Example
 *
 * \include snippet/core/type_traits/function_traits.cpp
 */
template <typename return_t, typename ...args_t>
struct function_traits<std::function<return_t(args_t...)>>
{
    //!\brief The number of arguments passed to the std::function target.
    static constexpr size_t argument_count = sizeof...(args_t);

    //!\brief The return type of the function target.
    using result_type = return_t;

    /*!\brief The argument type at the given `index`.
     * \tparam index The position of the argument to get the type for; must be smaller than `argument_count`.
     */
    template <size_t index>
    //!\cond
        requires index < argument_count
    //!\endcond
    using argument_type_at = pack_traits::at<index, args_t...>;
};

//!\cond
// Overload for all function types.
template <typename function_t>
    requires requires (function_t fn) { {std::function{fn}}; }
struct function_traits<function_t> : function_traits<decltype(std::function{std::declval<function_t>()})>
{};
//!\endcond
} // namespace seqan3

namespace seqan3::detail
{
// ----------------------------------------------------------------------------
// multi_invocable
// ----------------------------------------------------------------------------

/*!\brief A type that can conveniently inherit multiple invocables and acts as a union over them.
 * \tparam invocable_ts The types to inherit from.
 * \ingroup type_traits
 */
template <typename ...invocable_ts>
struct multi_invocable : invocable_ts...
{
    //!\brief Inherit the function call operators.
    using invocable_ts::operator()...;
};

//!\brief Deduction guides for seqan3::detail::multi_invocable.
template <typename ...invocable_ts>
multi_invocable(invocable_ts...) -> multi_invocable<invocable_ts...>;

} // namespace seqan3::detail
