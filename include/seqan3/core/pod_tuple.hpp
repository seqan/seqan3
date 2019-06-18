// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains seqan3::pod_tuple
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <meta/meta.hpp>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

//!\cond
#define SEQAN_NOT_POD "If you are not going to insert a POD type, use std::tuple instead."
//!\endcond

//!cond
template <typename ...types>
struct pod_tuple
{};
//!endcond

/*!\brief Behaves like std::tuple but is an aggregate [PODType](http://en.cppreference.com/w/cpp/concept/PODType).
 * \ingroup core
 * \implements seqan3::tuple_like_concept
 * \tparam type0    The first type (the first type).
 * \tparam ...types 0-n types (the remaining types of the values to be stored).
 *
 * This class behaves like std::tuple, but it is itself a POD type while std::tuple is not (even
 * if all contained types are POD). Since the only benefit of this class is that it stays POD it
 * actually enforces this on all types in the tuple (if you want to add non POD types, just use
 * std::tuple instead).
 *
 * It (only) supports [aggregate initialization](http://en.cppreference.com/w/cpp/language/aggregate_initialization),
 * i.e. you must use brace-initializiers and cannot
 * use paranthesis. You can use seqan3::get or std::get and also
 * [structured bindings](http://en.cppreference.com/w/cpp/language/declarations#Structured_binding_declaration)
 * to access the elements in the tuple.
 *
 * \snippet test/snippet/core/pod_tuple.cpp usage
 *
 */
template <typename type0, typename ...types>
struct pod_tuple<type0, types...>
{
    static_assert(std::is_pod_v<type0>, SEQAN_NOT_POD);
    //!\cond DEV
    //!\brief The first element as member.
    type0 _head;
    //!\brief The rest of the elements defined as a "recursive member".
    pod_tuple<types...> _tail;
    //!\endcond

    /*!\name Comparison operators
     * \{
     * \brief Lexicographically compares the values in the tuple.
     */

    //!\brief Test for equality.
    constexpr bool operator==(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) == std::tie(rhs._head, rhs._tail);
    }

    //!\brief Test for inequality.
    constexpr bool operator!=(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) != std::tie(rhs._head, rhs._tail);
    }

    //!\brief Test for smaller.
    constexpr bool operator<(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) < std::tie(rhs._head, rhs._tail);
    }

    //!\brief Test for larger.
    constexpr bool operator>(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) > std::tie(rhs._head, rhs._tail);
    }

    //!\brief Test for smaller or equal.
    constexpr bool operator<=(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) <= std::tie(rhs._head, rhs._tail);
    }

    //!\brief Test for larger or equal.
    constexpr bool operator>=(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) >= std::tie(rhs._head, rhs._tail);
    }
    //!\}
};

/*!\brief Recursion anchor for pod_tuple.
 * \ingroup core
 * \tparam type0 The value's type (every tuple must contain at least one type).
 */
template <typename type0>
struct pod_tuple<type0>
{
    static_assert(std::is_pod_v<type0>, SEQAN_NOT_POD);
    //!\cond DEV
    //!\brief The first element as member.
    type0 _head;
    //!\endcond

    /*!\name Comparison operators
     * \brief Lexicographically compares the values in the tuple.
     * \{
     */

    //!\brief Test for equality.
    constexpr bool operator==(pod_tuple const & rhs) const noexcept
    {
        return _head == rhs._head;
    }

    //!\brief Test for inequality.
    constexpr bool operator!=(pod_tuple const & rhs) const noexcept
    {
        return _head != rhs._head;
    }

    //!\brief Test for smaller.
    constexpr bool operator<(pod_tuple const & rhs) const noexcept
    {
        return _head < rhs._head;
    }

    //!\brief Test for larger.
    constexpr bool operator>(pod_tuple const & rhs) const noexcept
    {
        return _head > rhs._head;
    }

    //!\brief Test for smaller or equal.
    constexpr bool operator<=(pod_tuple const & rhs) const noexcept
    {
        return _head <= rhs._head;
    }

    //!\brief Test for larger or equal.
    constexpr bool operator>=(pod_tuple const & rhs) const noexcept
    {
        return _head >= rhs._head;
    }
    //!\}
};

#undef SEQAN_NOT_POD

//!\brief User defined deduction guide enables easy use.
//!\relates pod_tuple
template <typename ...types>
pod_tuple(types && ...) -> pod_tuple<types...>;

/*!\name Access an element of a pod_tuple by index
 * \{
 * \brief The same as [std::get](http://en.cppreference.com/w/cpp/utility/tuple/get) on an std::tuple.
 *
 * Note that these functions are available, both, in the seqan3 namespace and in namespace std.
 */
//!\relates seqan3::pod_tuple
template <std::size_t i, typename ...types>
constexpr auto & get(seqan3::pod_tuple<types...> & t) noexcept
    requires i < sizeof...(types)
{
    if constexpr (i == 0)
        return t._head;
    else
        return seqan3::get<i-1>(t._tail);
}

//!\relates seqan3::pod_tuple
template <std::size_t i, typename ...types>
constexpr auto const & get(seqan3::pod_tuple<types...> const & t) noexcept
    requires i < sizeof...(types)
{
    if constexpr (i == 0)
        return t._head;
    else
        return seqan3::get<i-1>(t._tail);
}

// extra overloads for temporaries required, because members of temporaries may only be returned as temporaries
//!\relates seqan3::pod_tuple
template <std::size_t i, typename ...types>
constexpr auto && get(seqan3::pod_tuple<types...> && t) noexcept
    requires i < sizeof...(types)
{
    if constexpr (i == 0)
        return std::move(t._head);
    else
        return seqan3::get<i-1>(std::move(t._tail));
}

//!\relates seqan3::pod_tuple
template <std::size_t i, typename ...types>
constexpr auto const && get(seqan3::pod_tuple<types...> const && t) noexcept
    requires i < sizeof...(types)
{
    if constexpr (i == 0)
        return std::move(t._head);
    else
        return seqan3::get<i-1>(std::move(t._tail));
}
//!\}

/*!\name Access an element of a pod_tuple by type
 * \brief The same as [std::get](http://en.cppreference.com/w/cpp/utility/tuple/get) on an std::tuple.
 *
 * Note that these functions are available, both, in the seqan3 namespace and in namespace std.
 * As is the case with std::tuple, this function is only defined if the type appears once
 * in the tuple, i.e. `std::get<int>(std::tuple<int, int>{1,2})` is not defined.
 * \{
 */

//!\relates seqan3::pod_tuple
template <typename type, typename ...types>
constexpr auto & get(seqan3::pod_tuple<types...> & t) noexcept
    requires meta::in<meta::list<types...>, type>::value &&
             (meta::find_index<meta::list<types...>, type>::value ==
              meta::reverse_find_index<meta::list<types...>, type>::value)
{
    return seqan3::get<meta::find_index<meta::list<types...>, type>::value>(t);
}

//!\relates seqan3::pod_tuple
template <typename type, typename ...types>
constexpr auto const & get(seqan3::pod_tuple<types...> const & t) noexcept
    requires meta::in<meta::list<types...>, type>::value &&
             (meta::find_index<meta::list<types...>, type>::value ==
              meta::reverse_find_index<meta::list<types...>, type>::value)
{
    return seqan3::get<meta::find_index<meta::list<types...>, type>::value>(t);
}

//!\relates seqan3::pod_tuple
template <typename type, typename ...types>
constexpr auto && get(seqan3::pod_tuple<types...> && t) noexcept
    requires meta::in<meta::list<types...>, type>::value &&
             (meta::find_index<meta::list<types...>, type>::value ==
              meta::reverse_find_index<meta::list<types...>, type>::value)
{
    return seqan3::get<meta::find_index<meta::list<types...>, type>::value>(std::move(t));
}

//!\relates seqan3::pod_tuple
template <typename type, typename ...types>
constexpr auto const && get(seqan3::pod_tuple<types...> const && t) noexcept
    requires meta::in<meta::list<types...>, type>::value &&
             (meta::find_index<meta::list<types...>, type>::value ==
              meta::reverse_find_index<meta::list<types...>, type>::value)
{
    return seqan3::get<meta::find_index<meta::list<types...>, type>::value>(std::move(t));
}
//!\}

} // namespace seqan3

namespace std
{

//!\cond
template <std::size_t i, typename ...types>
constexpr auto & get(seqan3::pod_tuple<types...> & t) noexcept
    requires i < sizeof...(types)
{
    return seqan3::get<i>(t);
}

template <std::size_t i, typename ...types>
constexpr auto const & get(seqan3::pod_tuple<types...> const & t) noexcept
    requires i < sizeof...(types)
{
    return seqan3::get<i>(t);
}

template <std::size_t i, typename ...types>
constexpr auto && get(seqan3::pod_tuple<types...> && t) noexcept
    requires i < sizeof...(types)
{
    return seqan3::get<i>(std::move(t));
}

template <std::size_t i, typename ...types>
constexpr auto const && get(seqan3::pod_tuple<types...> const && t) noexcept
    requires i < sizeof...(types)
{
    return seqan3::get<i>(std::move(t));
}

template <typename type, typename ...types>
constexpr auto & get(seqan3::pod_tuple<types...> & t) noexcept
    requires meta::in<meta::list<types...>, type>::value &&
             (meta::find_index<meta::list<types...>, type>::value ==
              meta::reverse_find_index<meta::list<types...>, type>::value)
{
    return seqan3::get<type>(t);
}

template <typename type, typename ...types>
constexpr auto const & get(seqan3::pod_tuple<types...> const & t) noexcept
    requires meta::in<meta::list<types...>, type>::value &&
             (meta::find_index<meta::list<types...>, type>::value ==
              meta::reverse_find_index<meta::list<types...>, type>::value)
{
    return seqan3::get<type>(t);
}

template <typename type, typename ...types>
constexpr auto && get(seqan3::pod_tuple<types...> && t) noexcept
    requires meta::in<meta::list<types...>, type>::value &&
             (meta::find_index<meta::list<types...>, type>::value ==
              meta::reverse_find_index<meta::list<types...>, type>::value)
{
    return seqan3::get<type>(std::move(t));
}

template <typename type, typename ...types>
constexpr auto const && get(seqan3::pod_tuple<types...> const && t) noexcept
    requires meta::in<meta::list<types...>, type>::value &&
             (meta::find_index<meta::list<types...>, type>::value ==
              meta::reverse_find_index<meta::list<types...>, type>::value)
{
    return seqan3::get<type>(std::move(t));
}
//!\endcond

//!\brief Obtains the type of the specified element.
//!\relates seqan3::pod_tuple
template <std::size_t i, template <typename...> typename t, typename ...types >
    requires i < sizeof...(types) &&
            std::is_base_of_v<seqan3::pod_tuple<types...>, t<types...>>
struct tuple_element<i, t<types...>>
{
    using type = meta::at_c<meta::list<types...>, i>;
};

//!\brief Provides access to the number of elements in a tuple as a compile-time constant expression.
//!\relates seqan3::pod_tuple
template <template <typename...> typename t, typename ...types >
    requires std::is_base_of_v<seqan3::pod_tuple<types...>, t<types...>>
struct tuple_size<t<types...>> :
    public std::integral_constant<std::size_t, sizeof...(types)>
{};

} // namespace std
