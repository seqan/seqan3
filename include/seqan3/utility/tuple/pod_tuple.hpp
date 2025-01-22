// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::pod_tuple
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <seqan3/utility/concept.hpp>
#include <seqan3/utility/type_pack/traits.hpp>

namespace seqan3
{

//!\cond
#define SEQAN_NOT_POD "If you are not going to insert a POD type, use std::tuple instead."
//!\endcond

template <typename... types>
struct pod_tuple
{};

/*!\brief Behaves like std::tuple but is an aggregate [PODType](https://en.cppreference.com/w/cpp/named_req/PODType).
 * \ingroup utility_tuple
 * \implements seqan3::tuple_like
 * \tparam type0    The first type (the first type).
 * \tparam ...types 0-n types (the remaining types of the values to be stored).
 *
 * This class behaves like std::tuple, but it is itself a POD type while std::tuple is not (even
 * if all contained types are POD). Since the only benefit of this class is that it stays POD it
 * actually enforces this on all types in the tuple (if you want to add non POD types, just use
 * std::tuple instead).
 *
 * You can use seqan3::get or std::get and also
 * [structured bindings](https://en.cppreference.com/w/cpp/language/declarations#Structured_binding_declaration)
 * to access the elements in the tuple.
 *
 * \include test/snippet/utility/tuple/pod_tuple.cpp
 *
 */
template <typename type0, typename... types>
struct pod_tuple<type0, types...>
{
    static_assert(std::is_standard_layout_v<type0> && seqan3::trivial<type0>, SEQAN_NOT_POD);
    //!\cond DEV
    //!\brief The first element as member.
    type0 _head;
    //!\brief The rest of the elements defined as a "recursive member".
    pod_tuple<types...> _tail;

    constexpr pod_tuple() noexcept = default;                              //!< Defaulted.
    constexpr pod_tuple(pod_tuple const &) noexcept = default;             //!< Defaulted.
    constexpr pod_tuple & operator=(pod_tuple const &) noexcept = default; //!< Defaulted.
    constexpr pod_tuple(pod_tuple &&) noexcept = default;                  //!< Defaulted.
    constexpr pod_tuple & operator=(pod_tuple &&) noexcept = default;      //!< Defaulted.
    constexpr ~pod_tuple() noexcept = default;                             //!< Defaulted.

    //!\brief Construct from arguments.
    constexpr pod_tuple(type0 v0, types... args) noexcept : _head{v0}, _tail{args...}
    {}
    //!\endcond

    /*!\name Comparison operators
     * \{
     */

    //!\brief Checks whether `*this` is equal to `rhs`.
    constexpr bool operator==(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) == std::tie(rhs._head, rhs._tail);
    }

    //!\brief Checks whether `*this` is not equal to `rhs`.
    constexpr bool operator!=(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) != std::tie(rhs._head, rhs._tail);
    }

    //!\brief Checks whether `*this` is less than `rhs`.
    constexpr bool operator<(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) < std::tie(rhs._head, rhs._tail);
    }

    //!\brief Checks whether `*this` is greater than `rhs`.
    constexpr bool operator>(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) > std::tie(rhs._head, rhs._tail);
    }

    //!\brief Checks whether `*this` is less than or equal to `rhs`.
    constexpr bool operator<=(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) <= std::tie(rhs._head, rhs._tail);
    }

    //!\brief Checks whether `*this` is greater than or equal to `rhs`.
    constexpr bool operator>=(pod_tuple const & rhs) const noexcept
    {
        return std::tie(_head, _tail) >= std::tie(rhs._head, rhs._tail);
    }
    //!\}
};

/*!\brief Recursion anchor for seqan3::pod_tuple.
 * \ingroup utility_tuple
 * \tparam type0 The value's type (every tuple must contain at least one type).
 */
template <typename type0>
struct pod_tuple<type0>
{
    static_assert(std::is_standard_layout_v<type0> && seqan3::trivial<type0>, SEQAN_NOT_POD);
    //!\cond DEV
    //!\brief The first element as member.
    type0 _head;
    //!\endcond

    /*!\name Comparison operators
     * \brief Lexicographically compares the values in the tuple.
     * \{
     */

    //!\brief Checks whether `*this` is equal to `rhs`.
    constexpr bool operator==(pod_tuple const & rhs) const noexcept
    {
        return _head == rhs._head;
    }

    //!\brief Checks whether `*this` is not equal to `rhs`.
    constexpr bool operator!=(pod_tuple const & rhs) const noexcept
    {
        return _head != rhs._head;
    }

    //!\brief Checks whether `*this` is less than `rhs`.
    constexpr bool operator<(pod_tuple const & rhs) const noexcept
    {
        return _head < rhs._head;
    }

    //!\brief Checks whether `*this` is greater than `rhs`.
    constexpr bool operator>(pod_tuple const & rhs) const noexcept
    {
        return _head > rhs._head;
    }

    //!\brief Checks whether `*this` is less than or equal to `rhs`.
    constexpr bool operator<=(pod_tuple const & rhs) const noexcept
    {
        return _head <= rhs._head;
    }

    //!\brief Checks whether `*this` is greater than or equal to `rhs`.
    constexpr bool operator>=(pod_tuple const & rhs) const noexcept
    {
        return _head >= rhs._head;
    }
    //!\}
};

#undef SEQAN_NOT_POD

//!\brief User defined deduction guide enables easy use.
//!\relates pod_tuple
template <typename... types>
pod_tuple(types &&...) -> pod_tuple<types...>;

/*!\name Access an element of a pod_tuple by index
 * \{
 * \brief The same as [std::get](https://en.cppreference.com/w/cpp/utility/tuple/get) on a std::tuple.
 *
 * Note that these functions are available, both, in the seqan3 namespace and in namespace std.
 */
//!\brief The same as [std::get](https://en.cppreference.com/w/cpp/utility/tuple/get) on a std::tuple.
//!\relates seqan3::pod_tuple
template <std::size_t i, typename... types>
constexpr auto & get(seqan3::pod_tuple<types...> & t) noexcept
    requires (i < sizeof...(types))
{
    if constexpr (i == 0)
        return t._head;
    else
        return seqan3::get<i - 1>(t._tail);
}

//!\brief The same as [std::get](https://en.cppreference.com/w/cpp/utility/tuple/get) on a std::tuple.
//!\relates seqan3::pod_tuple
template <std::size_t i, typename... types>
constexpr auto const & get(seqan3::pod_tuple<types...> const & t) noexcept
    requires (i < sizeof...(types))
{
    if constexpr (i == 0)
        return t._head;
    else
        return seqan3::get<i - 1>(t._tail);
}

// extra overloads for temporaries required, because members of temporaries may only be returned as temporaries
//!\brief The same as [std::get](https://en.cppreference.com/w/cpp/utility/tuple/get) on a std::tuple.
//!\relates seqan3::pod_tuple
template <std::size_t i, typename... types>
constexpr auto && get(seqan3::pod_tuple<types...> && t) noexcept
    requires (i < sizeof...(types))
{
    if constexpr (i == 0)
        return std::move(t._head);
    else
        return seqan3::get<i - 1>(std::move(t._tail));
}

//!\brief The same as [std::get](https://en.cppreference.com/w/cpp/utility/tuple/get) on a std::tuple.
//!\relates seqan3::pod_tuple
template <std::size_t i, typename... types>
constexpr auto const && get(seqan3::pod_tuple<types...> const && t) noexcept
    requires (i < sizeof...(types))
{
    if constexpr (i == 0)
        return std::move(t._head);
    else
        return seqan3::get<i - 1>(std::move(t._tail));
}
//!\}

/*!\name Access an element of a pod_tuple by type
 * \brief The same as [std::get](https://en.cppreference.com/w/cpp/utility/tuple/get) on a std::tuple.
 *
 * Note that these functions are available, both, in the seqan3 namespace and in namespace std.
 * As is the case with std::tuple, this function is only defined if the type appears once
 * in the tuple, i.e. `std::get<int>(std::tuple<int, int>{1,2})` is not defined.
 * \{
 */
//!\brief The same as [std::get](https://en.cppreference.com/w/cpp/utility/tuple/get) on a std::tuple.
//!\relates seqan3::pod_tuple
template <typename type, typename... arg_types>
constexpr auto & get(seqan3::pod_tuple<arg_types...> & t) noexcept
    requires (seqan3::pack_traits::count<type, arg_types...> == 1)
{
    return seqan3::get<seqan3::pack_traits::find<type, arg_types...>>(t);
}

//!\brief The same as [std::get](https://en.cppreference.com/w/cpp/utility/tuple/get) on a std::tuple.
//!\relates seqan3::pod_tuple
template <typename type, typename... arg_types>
constexpr auto const & get(seqan3::pod_tuple<arg_types...> const & t) noexcept
    requires (seqan3::pack_traits::count<type, arg_types...> == 1)
{
    return seqan3::get<seqan3::pack_traits::find<type, arg_types...>>(t);
}

//!\brief The same as [std::get](https://en.cppreference.com/w/cpp/utility/tuple/get) on a std::tuple.
//!\relates seqan3::pod_tuple
template <typename type, typename... arg_types>
constexpr auto && get(seqan3::pod_tuple<arg_types...> && t) noexcept
    requires (seqan3::pack_traits::count<type, arg_types...> == 1)
{
    return seqan3::get<seqan3::pack_traits::find<type, arg_types...>>(std::move(t));
}

//!\brief The same as [std::get](https://en.cppreference.com/w/cpp/utility/tuple/get) on a std::tuple.
//!\relates seqan3::pod_tuple
template <typename type, typename... arg_types>
constexpr auto const && get(seqan3::pod_tuple<arg_types...> const && t) noexcept
    requires (seqan3::pack_traits::count<type, arg_types...> == 1)
{
    return seqan3::get<seqan3::pack_traits::find<type, arg_types...>>(std::move(t));
}
//!\}

} // namespace seqan3

namespace std
{

//!\cond
template <std::size_t i, typename... types>
constexpr auto & get(seqan3::pod_tuple<types...> & t) noexcept
    requires (i < sizeof...(types))
{
    return seqan3::get<i>(t);
}

template <std::size_t i, typename... types>
constexpr auto const & get(seqan3::pod_tuple<types...> const & t) noexcept
    requires (i < sizeof...(types))
{
    return seqan3::get<i>(t);
}

template <std::size_t i, typename... types>
constexpr auto && get(seqan3::pod_tuple<types...> && t) noexcept
    requires (i < sizeof...(types))
{
    return seqan3::get<i>(std::move(t));
}

template <std::size_t i, typename... types>
constexpr auto const && get(seqan3::pod_tuple<types...> const && t) noexcept
    requires (i < sizeof...(types))
{
    return seqan3::get<i>(std::move(t));
}

template <typename type, typename... types>
constexpr auto & get(seqan3::pod_tuple<types...> & t) noexcept
    requires (seqan3::pack_traits::count<type, types...> == 1)
{
    return seqan3::get<type>(t);
}

template <typename type, typename... types>
constexpr auto const & get(seqan3::pod_tuple<types...> const & t) noexcept
    requires (seqan3::pack_traits::count<type, types...> == 1)
{
    return seqan3::get<type>(t);
}

template <typename type, typename... types>
constexpr auto && get(seqan3::pod_tuple<types...> && t) noexcept
    requires (seqan3::pack_traits::count<type, types...> == 1)
{
    return seqan3::get<type>(std::move(t));
}

template <typename type, typename... types>
constexpr auto const && get(seqan3::pod_tuple<types...> const && t) noexcept
    requires (seqan3::pack_traits::count<type, types...> == 1)
{
    return seqan3::get<type>(std::move(t));
}
//!\endcond

/*!\brief Obtains the type of the specified element.
 * \implements seqan3::transformation_trait
 * \relates seqan3::pod_tuple
 * \see [std::tuple_element](https://en.cppreference.com/w/cpp/utility/tuple/tuple_element)
 */
template <std::size_t i, template <typename...> typename t, typename... types>
    requires (i < sizeof...(types)) && std::is_base_of_v<seqan3::pod_tuple<types...>, t<types...>>
struct tuple_element<i, t<types...>>
{
    //!\brief Element type.
    using type = seqan3::pack_traits::at<i, types...>;
};

/*!\brief Provides access to the number of elements in a tuple as a compile-time constant expression.
 * \implements seqan3::unary_type_trait
 * \see std::tuple_size_v
 * \relates seqan3::pod_tuple
 */
template <template <typename...> typename t, typename... types>
    requires std::is_base_of_v<seqan3::pod_tuple<types...>, t<types...>>
struct tuple_size<t<types...>> : public std::integral_constant<std::size_t, sizeof...(types)>
{};

} // namespace std
