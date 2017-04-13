// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file core/pod_tuple.hpp
 * \ingroup core
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains seqan3::pod_tuple
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <meta/meta.hpp>

namespace seqan3
{

//!\cond
#define SEQAN_NOT_POD "If you are not going to insert a POD type, use std::tuple instead."
//!\endcond

/*!\brief Behaves like std::tuple but is an aggregate [PODType](http://en.cppreference.com/w/cpp/concept/PODType).
 * \ingroup core
 * \tparam type0 The first value's type (every tuple must contain at least one type).
 * \tparam ...types 0-n further types (the types of the other values).
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
 * ~~~~~~~~~~~~~~~{.cpp}
 *
 * pod_tuple<int, float> t{3, 4.7};
 * static_assert(std::is_pod_v<pod_tuple<int, float>>);
 *
 * // template parameters are automatically deduced:
 * pod_tuple t2{17, 3.7f, 19l};
 *
 * std::cout << std::get<0>(t2) << '\n'; // 17
 *
 * auto [ i, f, l ] = t2; // creates an int i with value 17, float f...
 *
 * ~~~~~~~~~~~~~~~
 *
 */
template <typename type0, typename ...types>
struct pod_tuple
{
    static_assert(std::is_pod_v<type0>, SEQAN_NOT_POD);
    //!\cond DEV
    //!\brief The first element as member.
    type0 _head;
    //!\brief The rest of the elements defined as a "recursive member".
    pod_tuple<types...> _tail;
    //!\endcond

    //!\name Comparison operators
    //!\{
    //!\brief Lexicographically compares the values in the tuple.
    constexpr bool operator==(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) == std::tie(rhs._head, rhs._tail);
    }

    constexpr bool operator!=(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) != std::tie(rhs._head, rhs._tail);
    }

    constexpr bool operator<(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) < std::tie(rhs._head, rhs._tail);
    }

    constexpr bool operator>(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) > std::tie(rhs._head, rhs._tail);
    }

    constexpr bool operator<=(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) <= std::tie(rhs._head, rhs._tail);
    }

    constexpr bool operator>=(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) >= std::tie(rhs._head, rhs._tail);
    }
    //!\}
};

template <typename type0>
struct pod_tuple<type0>
{
    static_assert(std::is_pod_v<type0>, SEQAN_NOT_POD);
    type0 _head;

    constexpr bool operator==(pod_tuple const & rhs) const noexcept
    {
        return _head == rhs._head;
    }

    constexpr bool operator!=(pod_tuple const & rhs) const noexcept
    {
        return _head != rhs._head;
    }

    constexpr bool operator<(pod_tuple const & rhs) const noexcept
    {
        return _head < rhs._head;
    }

    constexpr bool operator>(pod_tuple const & rhs) const noexcept
    {
        return _head > rhs._head;
    }

    constexpr bool operator<=(pod_tuple const & rhs) const noexcept
    {
        return _head <= rhs._head;
    }

    constexpr bool operator>=(pod_tuple const & rhs) const noexcept
    {
        return _head >= rhs._head;
    }
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
