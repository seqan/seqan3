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

/*!\brief Implementation of the unary type trait for seqan3::is_function.
 * \ingroup type_traits
 *
 * \tparam t The type to check if it is a valid function, lambda function, std::function or function pointer type.
 *
 * \details
 *
 * Provides the member constant `value` which is set to `true` if one of the following expressions evaluates to `true`:
 * * Evaluating the type in std::is_function returns `true`
 * * Creating a std::function object from an instance of `t` is not ill-formed and causes no substitution error.
 * Otherwise, `value` is `false`.
 */
template <typename t>
struct is_function
{
private:
    /*!\brief Helper function to check if `t` is a valid function type.
     * \tparam fn_t The alleged function type.
     *
     * \param[in] dummy Dummy variable used to enable overload resolution.
     *
     * \details
     *
     * Tries to instantiate a std::function object from an instance of `fn_t`. If the expression is ill-formed this
     * overload is removed by the compiler from the valid overloads through SFINAE and uses the overload with the
     * ellipse parameter which has always the lowest ranking for overload resolution.
     *
     * \returns `true` if creating an instance of std::function from `fn_t` is valid, otherwise `false`.
     */
    template <typename fn_t>
    constexpr static auto is_convertible_to_function(void * SEQAN3_DOXYGEN_ONLY(dummy))
    //!\cond
        -> decltype(static_cast<void>(std::function{std::declval<fn_t>()}), bool{})
    //!\endcond
    {
        return true;
    }

    //!\overload
    template <typename fn_t>
    constexpr static bool is_convertible_to_function(...)
    {
        return false;
    }

public:
    //!\brief Member constant to indicate whether `t` is a function type.
    static constexpr bool value = std::is_function_v<t> || is_convertible_to_function<t>(static_cast<void *>(nullptr));
};

} // namespace seqan3::detail

namespace seqan3
{
// ----------------------------------------------------------------------------
// is_function
// ----------------------------------------------------------------------------

/*!\brief Checks if a type is a function type, a lambda function type, a std::function type or a function pointer type.
 * \ingroup type_traits
 * \implements unary_type_trait
 *
 * \tparam t The type to check if it is a valid function, lambda function, std::function or function pointer type.
 *
 * \details
 *
 * This unary type trait extends the std::is_function type trait by also checking if the type is a
 * lambda function type, a std::function type or a function pointer type. Provides the member constant `value` which is
 * equal to `true`, if `t` is one of the aforementioned types. Otherwise, `value` is `false`.
 *
 * \see seqan3::is_function_v
 *
 * ### Example
 *
 * \include test/snippet/core/type_traits/is_function_trait.cpp
 */
template <typename t>
struct is_function : std::bool_constant<detail::is_function<t>::value>
{};

/*!\brief Helper variable template for seqan3::is_function unary transformation trait.
 * \ingroup type_traits
 *
 * \tparam t The type to check if it is a valid function, lambda function, std::function or function pointer type.
 */
template <typename t>
inline constexpr bool is_function_v = is_function<t>::value;

// ----------------------------------------------------------------------------
// function_traits
// ----------------------------------------------------------------------------

//!\cond
template <typename function_t>
struct function_traits;
//!\endcond

/*!\brief A traits class to provide a uniform interface to the properties of a std::function type.
 * \ingroup type_traits
 *
 * seqan3::function_traits is the type trait class that provides a uniform interface to the properties of
 * a std::function type.
 * This makes it possible to access the return type and the argument types of the stored target function.
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

} // namespace seqan3
