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
// Author: Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
// ============================================================================

#pragma once

#include <tuple>
#include <type_traits>

#include <meta/meta.hpp>

/*!\file core/pod_tuple.hpp
 * \ingroup core
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Contains seqan3::pod_tuple
 */

namespace seqan3
{

//!\cond
#define SEQAN_NOT_POD "If you are not going to insert a POD type, use std::tuple instead."
//!\endcond

/*!\brief Behaves like std::tuple but std::is_pod and std::is_aggregate.
 * \ingroup core
 * \tparam type0 The first value's type (every tuple must contain at least one type).
 * \tparam ...types 0-n further types (the types of the other values).
 *
 * This class behaves like std::tuple, but it is itself a POD type which std::tuple is not, even
 * if all contained types are POD. Since the only benefit of this class is that it stays POD it
 * actually enforces this on all types in the tuple (if you want to add non POD types, just use
 * std::tuple instead).
 *
 * It (only) supports aggregate initialization, i.e. you must use brace-initializiers and cannot
 * use paranthesis.
 *
 * ~~~~~~~~~~~~~~~{.cpp}
 *
 * pod_tuple<int, float> t{3, 4.7};
 * static_assert(std::is_pod_v<pod_tuple<int, float>>);
 *
 * // template parameters are automatically deduced:
 * pod_tuple t2{17, 3.7f, 19l};
 *
 * ~~~~~~~~~~~~~~~
 *
 * \todo implement type-based std::get, make_pod_tuple
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

    //!\name comparison operators
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

//!\brief User defined deduction guide enables easy use.
template <typename ...types>
pod_tuple(types && ...) -> pod_tuple<types...>;

#undef SEQAN_NOT_POD

} // namespace seqan3

namespace std
{

/*!\name Access an element of a pod_tuple
 * \{
 * \brief The same as std::get on an std::tuple
 */
//!\relates seqan3::pod_tuple
template <std::size_t i, typename ...types>
constexpr auto & get(seqan3::pod_tuple<types...> & t)
    requires i < sizeof...(types)
{
    if constexpr (i == 0)
        return t._head;
    else
        return std::get<i-1>(t._tail);
}

//!\relates seqan3::pod_tuple
template <std::size_t i, typename ...types>
constexpr auto const & get(seqan3::pod_tuple<types...> const & t) noexcept
    requires i < sizeof...(types)
{
    if constexpr (i == 0)
        return t._head;
    else
        return std::get<i-1>(t._tail);
}

// extra overloads for temporaries required, because members of temporaries may only be returned as temporaries
//!\relates seqan3::pod_tuple
template <std::size_t i, typename ...types>
constexpr auto && get(seqan3::pod_tuple<types...> && t)
    requires i < sizeof...(types)
{
    if constexpr (i == 0)
        return std::move(t._head);
    else
        return std::move(std::get<i-1>(t._tail));
}

//!\relates seqan3::pod_tuple
template <std::size_t i, typename ...types>
constexpr auto const && get(seqan3::pod_tuple<types...> const && t) noexcept
    requires i < sizeof...(types)
{
    if constexpr (i == 0)
        return std::move(t._head);
    else
        return std::move(std::get<i-1>(t._tail));
}
//!\}

//!\brief Obtains the type of the specified element.
//!\relates seqan3::pod_tuple
template <std::size_t i, typename ...types >
    requires i < sizeof...(types)
struct tuple_element<i, seqan3::pod_tuple<types...>>
{
    using type = meta::at_c<meta::list<types...>, i>;
};

//!\brief Provides access to the number of elements in a tuple as a compile-time constant expression.
//!\relates seqan3::pod_tuple
template <typename ...types>
struct tuple_size<seqan3::pod_tuple<types...>> :
    public std::integral_constant<std::size_t, sizeof...(types)>
{};

} // namespace std
